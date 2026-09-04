use cubecl::{cube, CubeType};
use cubecl::frontend::Float;
use cubecl::frontend::InverseSqrtExpand;

#[derive(CubeType, Clone, Copy)]
pub struct Vec3<F: Float> {
    pub(crate) x: F,
    pub(crate) y: F,
    pub(crate) z: F,
}

#[cube]
pub fn vec3_cross<F: Float>(
    a: &Vec3<F>,
    b: &Vec3<F>
) -> Vec3<F> {
    Vec3::<F> {
        x: a.y * b.z - a.z * b.y,
        y: a.z * b.x - a.x * b.z,
        z: a.x * b.y - a.y * b.x,
    }
}

#[cube]
pub fn vec3_scalar_multiply<F: Float>(
    scalar: F,
    v: &Vec3<F>
) -> Vec3<F> {
    Vec3::<F> {
        x: v.x * scalar,
        y: v.y * scalar,
        z: v.z * scalar
    }
}

#[cube]
pub fn vec3_add<F: Float>(
    a: &Vec3<F>,
    b: &Vec3<F>,
) -> Vec3<F> {
    Vec3::<F> {
        x: a.x + b.x,
        y: a.y + b.y,
        z: a.z + b.z,
    }
}

#[cube]
pub fn vec3_dot<F: Float>(
    a: &Vec3<F>,
    b: &Vec3<F>,
) -> F {
    a.x * b.x + a.y * b.y + a.z * b.z
}

#[cube]
pub fn vec3_norm<F: Float>(
    v: &Vec3<F>,
) -> Vec3<F> {
    let self_dot = vec3_dot(v, v);
    let inv_magnitude = self_dot.inverse_sqrt();
    vec3_scalar_multiply(inv_magnitude, v)
}