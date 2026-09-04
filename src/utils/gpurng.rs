use cubecl::cube;
use cubecl::prelude::*;

#[cube]
fn xorshift32(state: &mut u32) -> u32 {
    let mut x = *state;

    x ^= x << 13;
    x ^= x >> 17;
    x ^= x << 5;

    *state = x;
    x
}

#[cube]
fn hash(mut x: u32) -> u32 {
    x ^= x >> 16;
    x *= 0x7feb352d;
    x ^= x >> 15;
    x *= 0x446ca68b;
    x ^= x >> 16;
    x
}

#[cube]
pub fn rng_state(
    layer:u32,
    uc_x:u32,
    uc_y:u32,
    site:u32,
) -> u32 {
    let mut state = hash(layer);
    state ^= hash(uc_x);
    state ^= hash(uc_y);
    state ^= hash(site);
    state
}

#[cube]
fn rand_uniform(state: &mut u32) -> f32 {
    let x = xorshift32(state);

    // 2^-32, avoids exactly 0
    (x as f32 + 1.0f32) * 2.3283064365386963e-10f32
}

#[cube]
pub fn gaussian_pair(state: &mut u32) -> (f32, f32) {
    let u1 = rand_uniform(state);
    let u2 = rand_uniform(state);

    let r = (-2.0 * u1.ln()).sqrt();
    let theta = 6.283185307179586 * u2;

    (
        r * theta.cos(),
        r * theta.sin(),
    )
}