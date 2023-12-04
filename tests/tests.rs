use quickcheck::{Arbitrary, Gen};
use quickcheck_macros::quickcheck;
use vmmc::{
    position::Position,
    simbox::{self, SimBox},
    InputParams,
};

// x_idx * cells_per_axis[1] + y_idx
// indexes start at bottom left and scan bottom to top

// #[quickcheck]
// fn simbox_init(simbox: SimBox) -> bool {
//     true
// }

/*
dimensions = (40.000, 126.000)
cell_width = 1.2121212121212122, cell_height = 1.2 max_interaction_range = 1.188739079101486
cells_per_axis: [33, 105]

up left 105 104      // -1
up right 3570 3464   // - 106
bot left 0 0
bot right 3465 3359   // -106
*/

#[quickcheck]
fn simbox_4_corners(simbox: SimBox) -> bool {
    let epsilon = 0.001;
    let up_left = simbox.get_cell_id(Position::new([
        simbox.min_x() + epsilon,
        simbox.max_y() - epsilon,
    ]));
    let up_right = simbox.get_cell_id(Position::new([
        simbox.max_x() - epsilon,
        simbox.max_y() - epsilon,
    ]));
    let bottom_left = simbox.get_cell_id(Position::new([
        simbox.min_x() + epsilon,
        simbox.min_y() + epsilon,
    ]));
    let bottom_right = simbox.get_cell_id(Position::new([
        simbox.max_x() - epsilon,
        simbox.min_y() + epsilon,
    ]));

    let expected_up_left = simbox.cells_per_axis()[1] - 1;
    let expected_up_right = simbox.cells_per_axis()[0] * simbox.cells_per_axis()[1] - 1;
    let expected_bottom_left = 0;
    let expected_bottom_right = (simbox.cells_per_axis()[0] - 1) * simbox.cells_per_axis()[1];

    println!("up left {:?} {:?}", up_left, expected_up_left);
    println!("up right {:?} {:?}", up_right, expected_up_right);
    println!("bot left {:?} {:?}", bottom_left, expected_bottom_left);
    println!("bot right {:?} {:?}", bottom_right, expected_bottom_right);

    up_left == expected_up_left
        && up_right == expected_up_right
        && bottom_left == expected_bottom_left
        && bottom_right == expected_bottom_right
}

//fn f64_in_range()
