pub fn get_mag_neel(
    spins: &[f32],
    uc_dims:[usize;3],
    mask: &[u32],
) -> ([f64;3], [f64;3]) {
    // process as f64 for precision
    let mut total_mag:[f64;3] = [0.0; 3];
    let mut total_neel:[f64;3] = [0.0; 3];

    let atoms_per_layer = uc_dims[0]*uc_dims[1]*4;
    let layers = uc_dims[2]*6;
    let mut total_atoms = 0.0;

    // could be parallel, so far fast enough
    for layer in 0..layers {
        for atom in 0..atoms_per_layer {
            if mask[layer*atoms_per_layer + atom]==1 {
                total_atoms += 1.0;
                // odd layers +ve, even -ve for neel vector
                let neel_sign = ((layer%2)*2) as f64 -1.0;
                let index = 3*(layer*atoms_per_layer + atom);

                total_mag[0] += spins[index] as f64;
                total_mag[1] += spins[index+1] as f64;
                total_mag[2] += spins[index+2] as f64;

                total_neel[0] += spins[index] as f64 * neel_sign;
                total_neel[1] += spins[index+1] as f64 * neel_sign;
                total_neel[2] += spins[index+2] as f64 * neel_sign;
            }

        }
    }

    (
        [
            total_mag[0]/total_atoms,
            total_mag[1]/total_atoms,
            total_mag[2]/total_atoms,
        ],
        [
            total_neel[0]/total_atoms,
            total_neel[1]/total_atoms,
            total_neel[2]/total_atoms,
        ]
    )
}