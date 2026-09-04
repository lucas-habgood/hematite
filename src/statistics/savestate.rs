use std::fs::File;
use std::io;
use std::io::prelude::*;

pub fn save_state_to_file(
    output: &[f32],
    uc_dims: [usize; 3],
    filename: &str,
) -> io::Result<()> {
    let full_path = ["savestates/",filename].concat();
    let mut file = File::create(full_path)?;

    let uc_dims_bytes: Vec<u8> = uc_dims
        .iter().flat_map(|x| x.to_le_bytes()).collect();
    let spins_bytes: Vec<u8> = output
        .iter().flat_map(|x| x.to_le_bytes()).collect();

    file.write_all(&uc_dims_bytes)?;
    file.write_all(&spins_bytes)?;

    println!("Saved state to {}", filename);
    Ok(())
}

pub fn read_state_from_file(
    filename: &str,
) -> io::Result<([usize; 3], Vec<f32>)> {
    let full_path = ["savestates/",filename].concat();
    let mut file = File::open(full_path)?;
    let mut bytes = Vec::new();
    file.read_to_end(&mut bytes)?;

    // read the uc_dims
    let mut uc_dims:[usize;3] = [0;3];
    for i in 0..3 {
        let start = i * size_of::<usize>();
        let end = start + size_of::<usize>();

        let array: [u8; size_of::<usize>()] =
            bytes[start..end].try_into().unwrap();

        uc_dims[i] = usize::from_le_bytes(array);
    }
    // Read remaining bytes as f32s
    let spin_bytes = &bytes[3*size_of::<usize>()..];

    if spin_bytes.len() % size_of::<f32>() != 0 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "f32 data has an invalid length",
        ));
    }

    let output: Vec<f32> = spin_bytes
        .chunks_exact(4)
        .map(|chunk| {
            f32::from_le_bytes(chunk.try_into().unwrap())
        })
        .collect();

    println!("Loaded state of {} spins", output.len()/3);
    Ok((uc_dims, output))
}