use std::fs::OpenOptions;
use std::io::Write;

pub fn log_samples(
    filename: &str,
    samples:  &[f32],
    num_samples: usize,
    num_indices: usize,
) -> Result<(), Box<dyn std::error::Error>> {
    let mut file = OpenOptions::new().create(true).append(true).open(
        filename
    )?;
    for i in 0..num_samples {
        for j in 0..num_indices {
            let index = (i*num_indices + j)*3;
            write!(file, "{}, {}, {},    ",
                   samples[index], samples[index+1], samples[index+2])?;
        }
        writeln!(file)?;
    }

    Ok(())
}