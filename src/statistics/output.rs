use std::time::Instant;
use std::fs::OpenOptions;
use std::io::Write;
use crate::utils::config::config;

pub fn handle_output(
    iter:u32,
    start: &Instant,
    temperature:f32,
    charge: f32,
    mag:[f64;3],
    neel:[f64;3],
) -> Result<(), Box<dyn std::error::Error>> {
    let mut file = OpenOptions::new().create(true).append(true).open(
        config().output_file.clone()
    )?;
    let formatted_output = format!(
        "{}: {:.3?} \tT: {:.3} \tQ: {:.3} \tM: {:.5} {:.5} {:.5} \tN: {:.5} {:.5} {:.5}",
       iter, start.elapsed(), temperature, charge,
       mag[0], mag[1], mag[2], neel[0], neel[1], neel[2]);
    println!("{}", formatted_output);
    writeln!(file, "{}", formatted_output).expect("Failed to write to output file");
    Ok(())
}