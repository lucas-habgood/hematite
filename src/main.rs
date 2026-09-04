use cubecl::wgpu::WgpuRuntime;
use cubecl::cuda::CudaRuntime;
use rand::prelude::StdRng;
use rand::SeedableRng;
use crate::simulation::interactions::get_interactions;
use crate::simulation::launchgpu::launch_gpu;
use crate::simulation::lithomask::{allow_all_mask, mask_from_image};
use crate::statistics::savestate::read_state_from_file;
use crate::utils::config::config;
use crate::utils::cpuspinarray::{dimensions_to_unit_cells, rng_init_spin_array};
use crate::utils::createfiles::create_files_folders;

mod simulation {
    pub mod interactions;
    pub mod llg;
    pub mod field;
    pub mod accessgpuspins;
    pub mod launchgpu;
    pub mod unitcelldata;
    pub mod lithomask;
}
mod utils {
    pub mod vec3;
    pub mod cpuspinarray;
    pub mod constants;
    pub mod config;
    pub mod gpurng;
    pub mod createfiles;
}

mod visualisation {
    pub mod createplot;
    pub mod interpolategrid;
    pub mod plotgrid;
}

mod statistics {
    pub mod magvectors;
    pub mod output;
    pub mod logsamples;
    pub mod topological;
    pub mod savestate;
}

fn main() {
    println!("================================================================
██╗  ██╗███████╗███╗   ███╗ █████╗ ████████╗██╗████████╗███████╗
██║  ██║██╔════╝████╗ ████║██╔══██╗╚══██╔══╝██║╚══██╔══╝██╔════╝
███████║█████╗  ██╔████╔██║███████║   ██║   ██║   ██║   █████╗
██╔══██║██╔══╝  ██║╚██╔╝██║██╔══██║   ██║   ██║   ██║   ██╔══╝
██║  ██║███████╗██║ ╚═╝ ██║██║  ██║   ██║   ██║   ██║   ███████╗
╚═╝  ╚═╝╚══════╝╚═╝     ╚═╝╚═╝  ╚═╝   ╚═╝   ╚═╝   ╚═╝   ╚══════╝
Specialised hematite atomistic LLG spin evolution
Written in blazingly fast rust - to simulate rust
By Lucas Habgood
================================================================");
    create_files_folders().expect("Unable to create necessary files and folders");
    
    println!("Creating Spin Array");
    let mut rng = StdRng::seed_from_u64(config().rng_seed);


    let (uc_dims, spin_array, mask) = if config().reload_save {
        let (uc_dims, spin_array)=read_state_from_file(
            &*config().load_file).unwrap();
        let mask:Vec<u32> = mask_from_image(uc_dims, &*config().mask_file_name)
            .unwrap_or_else(|_| allow_all_mask(uc_dims));
        (uc_dims, spin_array, mask)
    }
    else {
        let uc_dims = dimensions_to_unit_cells(config().dimensions);
        let mask:Vec<u32> = mask_from_image(uc_dims, &*config().mask_file_name)
            .unwrap_or_else(|_| allow_all_mask(uc_dims));
        let (spin_array, num_mask_spins) = rng_init_spin_array(
            uc_dims, &mut rng, &mask
        );
        println!("Created {} Spins", num_mask_spins);
        (uc_dims, spin_array, mask)
    };

    println!("Creating Interactions");
    let interactions = get_interactions();
    println!("Preparing Kernel");
    if config().cuda {
        launch_gpu::<CudaRuntime>(uc_dims, interactions, spin_array, mask, &mut rng);
    }
    else {
        launch_gpu::<WgpuRuntime>(uc_dims, interactions, spin_array, mask, &mut rng);
    }
}
