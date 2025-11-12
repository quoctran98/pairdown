use std::fs::File;
use std::time::Instant;
use std::io::{BufReader, BufRead, BufWriter, Write, Error, ErrorKind, stdout};
use std::collections::HashSet;
use rand::{Rng, SeedableRng};
use rand_chacha::ChaCha8Rng;
use clap::Parser;

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
}

fn count_reads_and_validate_fastq(r1_fname: &str, r2_fname: &str) -> Result<usize, Error> {
    let mut r1_lines = BufReader::new(File::open(r1_fname)?).lines();
    let mut r2_lines = BufReader::new(File::open(r2_fname)?).lines();
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

fn choose_reads<R: Rng>(total: usize, kept: usize, rng: &mut R) -> Result<HashSet<usize>, Error> {
    if total < kept {
        return Err(Error::new(
            ErrorKind::InvalidInput,
            format!("Number of kept reads ({}) must be less than or equal to the total number of reads ({}) ", kept, total)
        ));
    }
    let kept_idx = rand::seq::index::sample(rng, total as usize, kept as usize);
    return Ok(HashSet::from_iter(kept_idx.into_iter()))
}

fn write_fastq(r1_old_fname: &str, r2_old_fname: &str, prefix: &str, kept_idx: HashSet<usize>) -> Result<(), Error> {
    // Set up the writer buffer
    let new_r1_fname = format!("{}{}", prefix, r1_old_fname);
    let new_r2_fname = format!("{}{}", prefix, r2_old_fname);
    let mut r1_writer = BufWriter::new(File::create(new_r1_fname)?);
    let mut r2_writer = BufWriter::new(File::create(new_r2_fname)?);

    // Manually loop to read FASTQs in chunks
    let mut r1_lines = BufReader::new(File::open(r1_old_fname)?).lines();
    let mut r2_lines = BufReader::new(File::open(r2_old_fname)?).lines();
    let mut curr_read_idx: usize = 0; // choose_reads() returns indices in the interval [0,readcount)
    loop {
        let r1_read: Vec<_> = r1_lines.by_ref().take(4).collect::<Result<_, _>>()?;
        let r2_read: Vec<_> = r2_lines.by_ref().take(4).collect::<Result<_, _>>()?;
        // We error checking here (that FASTQs end at the same time, etc.) because we did it before :)
        if r1_read.len() < 4 {
            break;
        }
        // If this read is in the kept index, add it
        if kept_idx.contains(&curr_read_idx) {
            for line in &r1_read {
                writeln!(r1_writer, "{}", line)?;
            }
            for line in &r2_read {
                writeln!(r2_writer, "{}", line)?;
            }
        }
        curr_read_idx += 1;
    }
    Ok(())
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
    let readcount = count_reads_and_validate_fastq(&args.read1, &args.read2).expect("Failed to read/validate FASTQ files.");
    
    // Choose reads to keep then write the new fastqs
    println!("\rFound {} reads. Writing {} reads to new FASTQ files...", &readcount, &args.keep);
    stdout().flush()?; 
    let kept_idx = choose_reads(readcount, args.keep, &mut rng);
    let _ = write_fastq(&args.read1, &args.read2, &output_prefix, kept_idx?);

    // Print out final stats :)
    let out_r1_fname = format!("{}{}", &output_prefix, &args.read1);
    let out_r2_fname = format!("{}{}", &output_prefix, &args.read2);
    println!("\rDownsampled {} from {} reads to {} and {} in {} seconds.", 
        args.keep, &readcount, out_r1_fname, out_r2_fname, start.elapsed().as_secs()
    );
    stdout().flush()?; 
    return Ok(())
}