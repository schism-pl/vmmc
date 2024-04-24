use rand::{rngs::SmallRng, SeedableRng};
use rand_distr::num_traits::Zero;
use vmmc::{
    morphology::Morphology,
    particle::Particle,
    position::{Orientation, Position},
    protocol::SynthesisProtocol,
    vmmc::Vmmc,
    vmmc_from_inputparams, InputParams,
};

fn empty_vmmc(shape: Morphology) -> Vmmc {
    let mut ip = InputParams::default();
    ip.seed = 1337;
    ip.sim_params.initial_particles = 0;
    ip.sim_params.shapes = vec![shape];
    ip.protocol = SynthesisProtocol::flat_protocol(0.0, 8.0, 20);

    let mut rng = SmallRng::seed_from_u64(ip.seed as u64);

    // Generate the simulator
    vmmc_from_inputparams(&ip, &mut rng)
}

fn two_particle_vmmc(
    shape: Morphology,
    p0: Position,
    or0: Orientation,
    p1: Position,
    or1: Orientation,
) -> Vmmc {
    let mut vmmc = empty_vmmc(shape);
    let part1 = Particle::new(0, p0, or0, 0);
    let part2 = Particle::new(1, p1, or1, 0);
    vmmc.simbox_mut().insert_particle(part1);
    vmmc.simbox_mut().insert_particle(part2);
    vmmc
}

// TODO: make sure pair potential works in both directions
// TODO: check vertical

#[test]
fn potential_test_3patch_1() {
    let shape = Morphology::regular_3patch(0.05);

    let p0 = Position::new([1.0, 1.0]);
    let or0 = Orientation::new([1.0, 0.0]);

    let p1 = Position::new([2.02, 1.0]);
    let or1 = Orientation::new([-1.0, 0.0]);

    let vmmc = two_particle_vmmc(shape, p0, or0, p1, or1);
    let energy01 = vmmc.compute_pair_energy(vmmc.particle(0), vmmc.particle(1));
    let energy10 = vmmc.compute_pair_energy(vmmc.particle(1), vmmc.particle(0));
    assert!(energy01.is_normal() || energy01.is_zero());
    assert!(energy10.is_normal() || energy10.is_zero());
    assert_eq!(energy10, -8.0);
    assert_eq!(energy01, -8.0);
}

#[test]
fn potential_test_3patch_2() {
    let shape = Morphology::regular_3patch(0.05);

    let p0 = Position::new([1.0, 1.0]);
    let or0 = Orientation::new([1.0, 0.0]);

    let p1 = Position::new([2.02, 1.0]);
    let or1 = Orientation::new([0.5, 3.0_f64.sqrt() / 2.0]);

    let vmmc = two_particle_vmmc(shape, p0, or0, p1, or1);
    let energy01 = vmmc.compute_pair_energy(vmmc.particle(0), vmmc.particle(1));
    let energy10 = vmmc.compute_pair_energy(vmmc.particle(1), vmmc.particle(0));
    assert!(energy01.is_normal() || energy01.is_zero());
    assert!(energy10.is_normal() || energy10.is_zero());
    assert_eq!(energy10, -8.0);
    assert_eq!(energy01, -8.0);
}

#[test]
fn potential_test_3patch_3() {
    let shape = Morphology::regular_3patch(0.05);

    let p0 = Position::new([1.0, 1.0]);
    let or0 = Orientation::new([1.0, 0.0]);

    let p1 = Position::new([2.02, 1.0]);
    let or1 = Orientation::new([0.5, -3.0_f64.sqrt() / 2.0]);
    println!("3patch_3: {or1}");

    let vmmc = two_particle_vmmc(shape, p0, or0, p1, or1);
    let energy01 = vmmc.compute_pair_energy(vmmc.particle(0), vmmc.particle(1));
    let energy10 = vmmc.compute_pair_energy(vmmc.particle(1), vmmc.particle(0));
    assert!(energy01.is_normal() || energy01.is_zero());
    assert!(energy10.is_normal() || energy10.is_zero());
    assert_eq!(energy10, -8.0);
    assert_eq!(energy01, -8.0);
}

#[test]
fn potential_test_4patch_1() {
    let shape = Morphology::regular_4patch(0.05);

    let p0 = Position::new([1.0, 1.0]);
    let or0 = Orientation::new([1.0, 0.0]);

    let p1 = Position::new([2.02, 1.0]);
    let or1 = Orientation::new([1.0, 0.0]);

    let vmmc = two_particle_vmmc(shape, p0, or0, p1, or1);
    let energy01 = vmmc.compute_pair_energy(vmmc.particle(0), vmmc.particle(1));
    let energy10 = vmmc.compute_pair_energy(vmmc.particle(1), vmmc.particle(0));
    assert!(energy01.is_normal() || energy01.is_zero());
    assert!(energy10.is_normal() || energy10.is_zero());
    assert_eq!(energy10, -8.0);
    assert_eq!(energy01, -8.0);
}

#[test]
fn potential_test_4patch_2() {
    let shape = Morphology::regular_4patch(0.05);

    let p0 = Position::new([1.0, 1.0]);
    let or0 = Orientation::new([1.0, 0.0]);

    let p1 = Position::new([2.02, 1.0]);
    let or1 = Orientation::new([0.0, 1.0]);

    let vmmc = two_particle_vmmc(shape, p0, or0, p1, or1);
    let energy01 = vmmc.compute_pair_energy(vmmc.particle(0), vmmc.particle(1));
    let energy10 = vmmc.compute_pair_energy(vmmc.particle(1), vmmc.particle(0));
    assert!(energy01.is_normal() || energy01.is_zero());
    assert!(energy10.is_normal() || energy10.is_zero());
    assert_eq!(energy10, -8.0);
    assert_eq!(energy01, -8.0);
}

#[test]
fn potential_test_4patch_3() {
    let shape = Morphology::regular_4patch(0.05);

    let p0 = Position::new([1.0, 1.0]);
    let or0 = Orientation::new([1.0, 0.0]);

    let p1 = Position::new([2.02, 1.0]);
    let or1 = Orientation::new([-1.0, 0.0]);

    let vmmc = two_particle_vmmc(shape, p0, or0, p1, or1);
    let energy01 = vmmc.compute_pair_energy(vmmc.particle(0), vmmc.particle(1));
    let energy10 = vmmc.compute_pair_energy(vmmc.particle(1), vmmc.particle(0));
    assert!(energy01.is_normal() || energy01.is_zero());
    assert!(energy10.is_normal() || energy10.is_zero());
    assert_eq!(energy10, -8.0);
    assert_eq!(energy01, -8.0);
}

#[test]
fn potential_test_4patch_4() {
    let shape = Morphology::regular_4patch(0.05);

    let p0 = Position::new([1.0, 1.0]);
    let or0 = Orientation::new([1.0, 0.0]);

    let p1 = Position::new([2.02, 1.0]);
    let or1 = Orientation::new([0.0, -1.0]);

    let vmmc = two_particle_vmmc(shape, p0, or0, p1, or1);
    let energy01 = vmmc.compute_pair_energy(vmmc.particle(0), vmmc.particle(1));
    let energy10 = vmmc.compute_pair_energy(vmmc.particle(1), vmmc.particle(0));
    assert!(energy01.is_normal() || energy01.is_zero());
    assert!(energy10.is_normal() || energy10.is_zero());
    assert_eq!(energy10, -8.0);
    assert_eq!(energy01, -8.0);
}

// #[quickcheck]
// fn simbox_4_corners(simbox: SimBox) -> bool {
//     let epsilon = 0.001;
//     let up_left = simbox.get_cell_id(Position::new([
//         simbox.min_x() + epsilon,
//         simbox.max_y() - epsilon,
//     ]));
//     let up_right = simbox.get_cell_id(Position::new([
//         simbox.max_x() - epsilon,
//         simbox.max_y() - epsilon,
//     ]));
//     let bottom_left = simbox.get_cell_id(Position::new([
//         simbox.min_x() + epsilon,
//         simbox.min_y() + epsilon,
//     ]));
//     let bottom_right = simbox.get_cell_id(Position::new([
//         simbox.max_x() - epsilon,
//         simbox.min_y() + epsilon,
//     ]));

//     let expected_up_left = simbox.cells_per_axis()[1] - 1;
//     let expected_up_right = simbox.cells_per_axis()[0] * simbox.cells_per_axis()[1] - 1;
//     let expected_bottom_left = 0;
//     let expected_bottom_right = (simbox.cells_per_axis()[0] - 1) * simbox.cells_per_axis()[1];

//     println!("up left {:?} {:?}", up_left, expected_up_left);
//     println!("up right {:?} {:?}", up_right, expected_up_right);
//     println!("bot left {:?} {:?}", bottom_left, expected_bottom_left);
//     println!("bot right {:?} {:?}", bottom_right, expected_bottom_right);

//     up_left == expected_up_left
//         && up_right == expected_up_right
//         && bottom_left == expected_bottom_left
//         && bottom_right == expected_bottom_right
// }

//fn f64_in_range()
