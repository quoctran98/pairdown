use std::fs::File;
use std::time::Instant;
use std::io::{BufReader, BufRead, BufWriter, Write, Error, ErrorKind, stdout};
use rand::{Rng, SeedableRng};
use rand_chacha::ChaCha8Rng;
use clap::Parser;

// const BUFCAPACITY: usize = 4 * 1024 * 1024; // 4 MB
// const BUFCAPACITY: usize = 8 * 1024; // 8 KB
const BUFCAPACITY: usize = 32 * 1024; // 32 KB


/// Downsample paired-end FASTQ files
#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// Path to Read 1 FASTQ file
    #[arg(index = 1)]
    read1: String,
    /// Path to Read 2 FASTQ file
    #[arg(index = 2)]
    read2: String,
    /// Number of reads to keep
    #[arg(index = 3)]
    keep: usize,
    /// Output prefix (default: downsampled_<reads_kept>_)
    #[arg(short = 'o', long)]
    prefix: Option<String>, // Option<> allows for this to be skipped 
    /// Random seed (default: 17)
    #[arg(short = 's', long, default_value_t = 17)]
    seed: u64,
    /// Validate FASTQs (much slower, but ensures FASTQs are paired)
    #[arg(short = 'v', long, default_value_t = false)]
    validate: bool,
}

// If we do trust that the two FASTQs are valid, we can just count a single one (faster but less safe)
// This version reuses the same read buffer vector to hold the current bytes
fn count_reads_fast(fname: &str) -> Result<usize, Error> {
    let mut reader = BufReader::with_capacity(BUFCAPACITY, File::open(fname)?);
    let mut buf = Vec::with_capacity(256); // We don't care about anything in here.
    let mut total_lines: usize = 0;
    while reader.read_until(b'\n', &mut buf)? != 0 {
        total_lines += 1;
        buf.clear();
    }
    return Ok(total_lines/4) // What happens if it's not divisible by four?
}

// If we don't trust that the two FASTQs are valid (same reads, lines divisible by 4)
fn count_reads_and_validate_fastq(r1_fname: &str, r2_fname: &str) -> Result<usize, Error> {
    let mut r1_lines = BufReader::with_capacity(BUFCAPACITY, File::open(r1_fname)?).lines();
    let mut r2_lines = BufReader::with_capacity(BUFCAPACITY, File::open(r2_fname)?).lines();
    let mut readcount: usize = 0;
    loop {
        // Read four lines at at time
        let r1_read: Vec<_> = r1_lines.by_ref().take(4).collect::<Result<_, _>>()?;
        let r2_read: Vec<_> = r2_lines.by_ref().take(4).collect::<Result<_, _>>()?;
        // Throw errors for invalid paired read FASTQs
        if r1_read.len() != r2_read.len() {
            return Err(Error::new(
                ErrorKind::InvalidData,
                "FASTQ files have different number of reads",
            ));
        }
        if r1_read.len() % 4 != 0 {
            return Err(Error::new(
                ErrorKind::InvalidData,
                format!("The number of lines in {} is not divisible by 4.", r1_fname)
            ));
        }
        if r2_read.len() % 4 != 0 {
            return Err(Error::new(
                ErrorKind::InvalidData,
                format!("The number of lines in {} is not divisible by 4.", r2_fname)
            ));
        }
        // Stop if we have no more lines; continue otherwise (error handling should have covered all other cases)
        if r1_read.len() == 0 {
            break;
        }
        readcount += 1;
    }
    return Ok(readcount)
}

fn choose_reads<R: Rng>(total: usize, kept: usize, rng: &mut R) -> Result<Vec<bool>, Error> {
    if total < kept {
        return Err(Error::new(
            ErrorKind::InvalidInput,
            format!("Number of kept reads ({}) must be less than or equal to the total number of reads ({}) ", kept, total)
        ));
    }
    let mut keep_vec = vec![false; total];
    for i in rand::seq::index::sample(rng, total as usize, kept as usize).iter() {
        keep_vec[i] = true;
    }
    return Ok(keep_vec)
}

fn write_fastq(r1_old_fname: &str, r2_old_fname: &str, prefix: &str, keep_vec: Vec<bool>) -> Result<(), Error> {
    // Set up the writer buffer
    let mut writer_r1 = BufWriter::with_capacity(BUFCAPACITY, File::create(format!("{}{}", prefix, r1_old_fname))?);
    let mut writer_r2 = BufWriter::with_capacity(BUFCAPACITY, File::create(format!("{}{}", prefix, r2_old_fname))?);

    // Read raw bytes from the FASTQ (without casting to Strings) and write the raw bytes to the new FASTQ
    let mut reader_r1 = BufReader::with_capacity(BUFCAPACITY, File::open(r1_old_fname)?);
    let mut reader_r2 = BufReader::with_capacity(BUFCAPACITY, File::open(r2_old_fname)?);
    let mut buf_r1 = Vec::with_capacity(512); // I should do the math to figure out max capacity (for 150 cycles or something)
    let mut buf_r2 = Vec::with_capacity(512);

    let mut curr_read_idx: usize = 0; // choose_reads() returns indices in the interval [0,readcount)
    'outer: loop {
        for _ in 0..4 { // Read until we've hit four newlines. Only break the big loop if we reach EOF
            if reader_r1.read_until(b'\n', &mut buf_r1)? == 0 {
                break 'outer;
            }
            if reader_r2.read_until(b'\n', &mut buf_r2)? == 0 {
                break 'outer;
            }
        }
        if keep_vec[curr_read_idx] { // Write the whole buffer if we're keeping this line.
            writer_r1.write_all(&buf_r1)?;
            writer_r2.write_all(&buf_r2)?;
        }
        buf_r1.clear();
        buf_r2.clear();
        curr_read_idx += 1;
    }
    return Ok(())
}

fn main() -> Result<(), Error> {

    // Keep track of time for the fun of it, I guess...
    let start = Instant::now();

    // Parse args and set up seeds
    let args = Args::parse();
    let output_prefix = args.prefix.unwrap_or_else(|| format!("downsampled_{}_", args.keep));
    let mut rng = ChaCha8Rng::seed_from_u64(args.seed);

    // Count reads in both FASTQs (this will also kind of validate the paired end read files)
    print!("Counting reads in {} and {} ...", &args.read1, &args.read2);
    stdout().flush()?; 
    // Choose the function dpending on if fast mode is here
    let readcount: usize;
    if args.validate {
        readcount = count_reads_and_validate_fastq(&args.read1, &args.read2).expect("Failed to read or validate FASTQ files.");
    } else {
        readcount = count_reads_fast(&args.read1).expect("Failed to count reads in FASTQ file (no validation).");
    }
    
    // Choose reads to keep then write the new fastqs
    println!("\rFound {} reads. Writing {} reads to new FASTQ files...", &readcount, &args.keep);
    stdout().flush()?; 
    let keep_vec = choose_reads(readcount, args.keep, &mut rng).expect("Failed to randomly sample read indicies.");
    let _ = write_fastq(&args.read1, &args.read2, &output_prefix, keep_vec);

    // Print out final stats :)
    let out_r1_fname = format!("{}{}", &output_prefix, &args.read1);
    let out_r2_fname = format!("{}{}", &output_prefix, &args.read2);
    println!("\rDownsampled {} from {} reads to {} and {} in {} seconds.", 
        args.keep, &readcount, out_r1_fname, out_r2_fname, start.elapsed().as_secs()
    );
    stdout().flush()?; 
    return Ok(())
}