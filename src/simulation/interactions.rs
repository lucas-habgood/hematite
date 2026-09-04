use cubecl::{CubeLaunch, CubeType};
use cubecl::prelude::*;
use cubecl::server::Handle;
use crate::simulation::unitcelldata::{UNIT_CELL_DIMENSIONS, UNIT_CELL_LAYERS_COORDINATES};
use crate::utils::config::config;

fn norm(v: [f32; 3]) -> f32 {
    (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()
}

fn rounded(x: f32) -> f32 {
    (x*100.0).round()/100.0
}
fn distance(u: [f32; 3], v:[f32; 3]) -> f32 {
    let difference: [f32; 3] = [u[0] - v[0], u[1] - v[1], u[2] - v[2]];
    rounded(norm(difference))
}

fn offset_vec(u: [f32; 3], v: [f32; 3]) -> [f32; 3] {
    [u[0] + v[0], u[1] + v[1], u[2] + v[2]]
}

pub struct Interaction {
    pub layer_offset: i32,
    pub uc_x_offset: i32,
    pub uc_y_offset: i32,
    pub pair_site: u32,
    pub exchange_tensor: ExchangeTensor
}

pub struct ExchangeTensor {
    pub jxx : f32,
    pub jxy : f32,
    pub jxz : f32,
    pub jyx : f32,
    pub jyy : f32,
    pub jyz : f32,
    pub jzx : f32,
    pub jzy : f32,
    pub jzz : f32,
}

fn zeros_interaction_tensor() -> ExchangeTensor {
    let tensor = ExchangeTensor {
        jxx:0.,
        jxy:0.,
        jxz:0.,
        jyx:0.,
        jyy:0.,
        jyz:0.,
        jzx:0.,
        jzy:0.,
        jzz:0.,
    };
    tensor
}

// symmetric exchange energies of 1st, 2nd, 3rd, 4th nearest neighbours
const SYMMETRIC_EXCHANGE_ENERGIES: [f32;4] = [1.08410544e-21, 1.8608832e-22, -7.15911372e-21,-4.42453176e-21];

fn apply_symmetric_exchange(
    tensor: &mut ExchangeTensor,
    energy: f32,
) {
    tensor.jxx=energy;
    tensor.jyy=energy;
    tensor.jzz=energy;
}

fn apply_dmi(
    tensor: &mut ExchangeTensor,
    dmi_vector: [f32; 3],
) {
    tensor.jxy = dmi_vector[2];
    tensor.jyx = -dmi_vector[2];
    tensor.jxz = -dmi_vector[1];
    tensor.jzx = dmi_vector[1];
    tensor.jyz = dmi_vector[0];
    tensor.jzy = -dmi_vector[0];
}

fn apply_interfacial_dmi(
    tensor: &mut ExchangeTensor,
    energy: f32,
    bond_normal: [f32;3],
    interface_normal: [f32;3],
) {
    let cross_product = [
        bond_normal[1]*interface_normal[2] - bond_normal[2]*interface_normal[1],
        bond_normal[2]*interface_normal[0] - bond_normal[0]*interface_normal[2],
        bond_normal[0]*interface_normal[1] - bond_normal[1]*interface_normal[0]
    ];
    let dmi_vector = [
        cross_product[0]*energy, cross_product[1]*energy, cross_product[2]*energy
    ];
    apply_dmi(tensor, dmi_vector);
}

pub fn get_interactions() -> PackedInteractions {
    let mut uc_neighbour_offsets: Vec<[i32;3]> = vec![];
    for i in -1i32..2 {
        for j in -1i32..2 {
            for k in -1i32..2 {
                uc_neighbour_offsets.push([i,j,k]);
            }
        }
    }

    let mut distances: Vec<f32> = vec![];
    let mut interactions: [[Vec<Interaction>;4];6] = Default::default();

    for layer in 0..6 { for site in 0..4 {
        let coord: [f32;3] = UNIT_CELL_LAYERS_COORDINATES[layer][site];

        for uc_offset in &uc_neighbour_offsets {
            let offset: [f32;3] = [
                (uc_offset[0] as f32)*UNIT_CELL_DIMENSIONS[0],
                (uc_offset[1] as f32)*UNIT_CELL_DIMENSIONS[1],
                (uc_offset[2] as f32)*UNIT_CELL_DIMENSIONS[2]
            ];

            for neighbour_layer in 0..6 { for neighbour_site in 0..4 {
                let neighbour_coord = offset_vec(
                    UNIT_CELL_LAYERS_COORDINATES[neighbour_layer][neighbour_site],
                    offset);
                let distance = distance(coord, neighbour_coord);
                if distance==0.0 {continue;}
                distances.push(distance);
                let bond_vector = [
                    neighbour_coord[0]-coord[0],
                    neighbour_coord[1]-coord[1],
                    neighbour_coord[2]-coord[2]
                ];
                let bond_length = (
                    bond_vector[0]*bond_vector[0]
                    + bond_vector[1]*bond_vector[1]
                    + bond_vector[2]*bond_vector[2]
                ).sqrt();
                let bond_normal = [
                    bond_vector[0]/bond_length,
                    bond_vector[1]/bond_length,
                    bond_vector[2]/bond_length
                ];

                let true_neighbour_layer = (neighbour_layer as i32)+uc_offset[2]*6;

                let mut interaction:Interaction = Interaction{
                    layer_offset: true_neighbour_layer - (layer as i32),
                    uc_x_offset: uc_offset[0],
                    uc_y_offset: uc_offset[1],
                    pair_site: neighbour_site as u32,
                    exchange_tensor: zeros_interaction_tensor()
                };

                let interfacial_dmi_energy = config().interfacial_dmi_energy;
                let effective_canting_dmi_energy = config().effective_canting_dmi_energy;
                let interface_normal = config().interface_normal;

                if distance==2.9 {
                    apply_symmetric_exchange(
                        &mut interaction.exchange_tensor, SYMMETRIC_EXCHANGE_ENERGIES[0]);
                    // for effective dmi, need to sum to overall z component considering spin
                    // so from even->odd has different sign odd->even
                    let sign = ((layer%2)*2) as f32 -1.0;
                    let dmi_vector = [0.,0.,effective_canting_dmi_energy*sign];
                    apply_dmi(&mut interaction.exchange_tensor, dmi_vector);
                    interactions[layer][site].push(interaction);
                }
                else if distance==3.42 {
                    apply_symmetric_exchange(
                        &mut interaction.exchange_tensor, SYMMETRIC_EXCHANGE_ENERGIES[2]);
                    apply_interfacial_dmi(
                        &mut interaction.exchange_tensor,
                        interfacial_dmi_energy,
                        bond_normal,
                        interface_normal
                    );
                    interactions[layer][site].push(interaction);
                }
                else if distance==3.75 {
                    apply_symmetric_exchange(
                        &mut interaction.exchange_tensor, SYMMETRIC_EXCHANGE_ENERGIES[3]);
                    apply_interfacial_dmi(
                        &mut interaction.exchange_tensor,
                        interfacial_dmi_energy,
                        bond_normal,
                        interface_normal
                    );
                    interactions[layer][site].push(interaction);
                }
                else if distance==3.0 && config().intralayer {
                    apply_symmetric_exchange(
                        &mut interaction.exchange_tensor, SYMMETRIC_EXCHANGE_ENERGIES[1]);
                    // no iDMI applied here - would give different iDMI/area when using or not
                    // bond is weak, maybe iDMI would be also
                    interactions[layer][site].push(interaction);
                }
            }}
        }
    }}

    distances.sort_by(f32::total_cmp);
    distances.dedup_by(|a, b| a.total_cmp(b).is_eq());

    pack_to_arrays(interactions)
}

pub struct PackedInteractions {
    pub layer_offset: Vec<i32>,
    pub uc_x_offset: Vec<i32>,
    pub uc_y_offset: Vec<i32>,
    pub pair_site: Vec<u32>,
    pub jxx: Vec<f32>,
    pub jxy: Vec<f32>,
    pub jxz: Vec<f32>,
    pub jyx: Vec<f32>,
    pub jyy: Vec<f32>,
    pub jyz: Vec<f32>,
    pub jzx: Vec<f32>,
    pub jzy: Vec<f32>,
    pub jzz: Vec<f32>,
}

pub fn pack_to_arrays(interactions:[[Vec<Interaction>; 4]; 6]) -> PackedInteractions {
    let mut packed = PackedInteractions {
        layer_offset: vec![],
        uc_x_offset: vec![],
        uc_y_offset: vec![],
        pair_site: vec![],
        jxx: vec![],
        jxy: vec![],
        jxz: vec![],
        jyx: vec![],
        jyy: vec![],
        jyz: vec![],
        jzx: vec![],
        jzy: vec![],
        jzz: vec![],
    };
    for layer in 0..6 {
        for site in 0..4 {
            let ints_vec = &interactions[layer][site];
            for pairing in 0..interactions[0][0].len() {
                let int = &ints_vec[pairing];
                packed.layer_offset.push(int.layer_offset);
                packed.uc_x_offset.push(int.uc_x_offset);
                packed.uc_y_offset.push(int.uc_y_offset);
                packed.pair_site.push(int.pair_site);
                packed.jxx.push(int.exchange_tensor.jxx);
                packed.jxy.push(int.exchange_tensor.jxy);
                packed.jxz.push(int.exchange_tensor.jxz);
                packed.jyx.push(int.exchange_tensor.jyx);
                packed.jyy.push(int.exchange_tensor.jyy);
                packed.jyz.push(int.exchange_tensor.jyz);
                packed.jzx.push(int.exchange_tensor.jzx);
                packed.jzy.push(int.exchange_tensor.jzy);
                packed.jzz.push(int.exchange_tensor.jzz);
            }
        }
    }

    packed
}


#[derive(CubeType, CubeLaunch)]
pub struct GPUInteractions {
    pub layer_offset: Array<i32>,
    pub uc_x_offset: Array<i32>,
    pub uc_y_offset: Array<i32>,
    pub pair_site: Array<u32>,
    pub jxx: Array<f32>,
    pub jxy: Array<f32>,
    pub jxz: Array<f32>,
    pub jyx: Array<f32>,
    pub jyy: Array<f32>,
    pub jyz: Array<f32>,
    pub jzx: Array<f32>,
    pub jzy: Array<f32>,
    pub jzz: Array<f32>,
}

pub struct GPUInteractionsStorage {
    pub layer_offset: Handle,
    pub uc_x_offset: Handle,
    pub uc_y_offset: Handle,
    pub pair_site: Handle,
    pub jxx: Handle,
    pub jxy: Handle,
    pub jxz: Handle,
    pub jyx: Handle,
    pub jyy: Handle,
    pub jyz: Handle,
    pub jzx: Handle,
    pub jzy: Handle,
    pub jzz: Handle,
}

pub fn upload_interactions<R: Runtime>(
    packed: &PackedInteractions,
    client: &ComputeClient<R>,
) -> GPUInteractionsStorage {
    GPUInteractionsStorage {
        layer_offset: client.create_from_slice(
            i32::as_bytes(&packed.layer_offset)
        ),
        uc_x_offset: client.create_from_slice(
            i32::as_bytes(&packed.uc_x_offset)
        ),
        uc_y_offset: client.create_from_slice(
            i32::as_bytes(&packed.uc_y_offset)
        ),
        pair_site: client.create_from_slice(
            u32::as_bytes(&packed.pair_site)
        ),
        jxx: client.create_from_slice(
            f32::as_bytes(&packed.jxx)
        ),
        jxy: client.create_from_slice(
            f32::as_bytes(&packed.jxy)
        ),
        jxz: client.create_from_slice(
            f32::as_bytes(&packed.jxz)
        ),
        jyx: client.create_from_slice(
            f32::as_bytes(&packed.jyx)
        ),
        jyy: client.create_from_slice(
            f32::as_bytes(&packed.jyy)
        ),
        jyz: client.create_from_slice(
            f32::as_bytes(&packed.jyz)
        ),
        jzx: client.create_from_slice(
            f32::as_bytes(&packed.jzx)
        ),
        jzy: client.create_from_slice(
            f32::as_bytes(&packed.jzy)
        ),
        jzz: client.create_from_slice(
            f32::as_bytes(&packed.jzz)
        ),
    }
}

pub unsafe fn interaction_launch<R: Runtime>(
    data: &GPUInteractionsStorage,
    num: usize,
) -> GPUInteractionsLaunch<R> { unsafe {
    GPUInteractionsLaunch::new(
        ArrayArg::from_raw_parts(data.layer_offset.clone(), 24*num),
        ArrayArg::from_raw_parts(data.uc_x_offset.clone(), 24*num),
        ArrayArg::from_raw_parts(data.uc_y_offset.clone(), 24*num),
        ArrayArg::from_raw_parts(data.pair_site.clone(), 24*num),
        ArrayArg::from_raw_parts(data.jxx.clone(), 24*num),
        ArrayArg::from_raw_parts(data.jxy.clone(), 24*num),
        ArrayArg::from_raw_parts(data.jxz.clone(), 24*num),
        ArrayArg::from_raw_parts(data.jyx.clone(), 24*num),
        ArrayArg::from_raw_parts(data.jyy.clone(), 24*num),
        ArrayArg::from_raw_parts(data.jyz.clone(), 24*num),
        ArrayArg::from_raw_parts(data.jzx.clone(), 24*num),
        ArrayArg::from_raw_parts(data.jzy.clone(), 24*num),
        ArrayArg::from_raw_parts(data.jzz.clone(), 24*num),
    )
}}