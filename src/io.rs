use std::{
    f32::consts::PI,
    fs::{self, File},
    io::{BufRead, BufReader, Write},
};

// use raqote::{PathBuilder, DrawTarget, Source, SolidSource, StrokeStyle, DrawOptions};
use raqote::*;
use crate::types::Num;
use crate::{cli::VmmcConfig, particle::IsParticle, types::DimVec, vmmc::Vmmc};

pub struct XYZWriter {
    file: File,
}

// pub struct XYZReader {
//     rdr: BufReader<File>,
// }

// impl XYZReader {
//     pub fn new(p: &str) -> Self {
//         let file = File::open(p).unwrap();
//         let rdr = BufReader::new(file);
//         Self { rdr }
//     }
// }

pub fn read_xyz_snapshot(path: &str) -> (Vec<DimVec>, Vec<DimVec>) {
    let mut positions = Vec::new();
    let mut orientations = Vec::new();

    let file = File::open(path).unwrap();
    let rdr = BufReader::new(file);

    for line in rdr.lines() {
        let line = line.unwrap();
        let parts: Vec<&str> = line.split_whitespace().collect();
        let pos_x = parts[0].parse::<Num>().unwrap();
        let pos_y = parts[1].parse::<Num>().unwrap();
        let or_x = parts[2].parse::<Num>().unwrap();
        let or_y = parts[3].parse::<Num>().unwrap();
        positions.push(DimVec::new([pos_x, pos_y]));
        orientations.push(DimVec::new([or_x, or_y]));
    }

    (positions, orientations)
}

impl XYZWriter {
    pub fn new(p: &str) -> Self {
        let file = File::create(p).unwrap();
        Self { file }
    }

    pub fn write_xyz_frame(&mut self, vmmc: &Vmmc) {
        writeln!(self.file, "{:?}\n", vmmc.particles().num_particles()).unwrap();
        for p in vmmc.particles().iter() {
            writeln!(self.file, "0 {:?} {:?} 0", p.pos().x(), p.pos().y()).unwrap();
        }
    }
}

pub fn write_tcl(vmmc: &Vmmc, p: &str) {
    let mut file = File::create(p).unwrap();
    write!(
        file,
        "light 0 on
    light 1 on
    light 2 off
    light 3 off
    axes location off
    stage location off
    display projection orthographic
    mol modstyle 0 0 VDW 1 30
    set sel [atomselect top \"name X\"]
    atomselect0 set radius 0.5
    color Name X blue
    display depthcue off
    set minx {:?}
    set maxx {:?}
    set miny {:?}
    set maxy {:?}
    set minz 0
    set maxz 0
    draw materials off
    draw color white
    draw line \"$minx $miny $minz\" \"$maxx $miny $minz\" 
    draw line \"$minx $miny $minz\" \"$minx $maxy $minz\" 
    draw line \"$minx $maxy $minz\" \"$maxx $maxy $minz\" 
    draw line \"$maxx $miny $minz\" \"$maxx $maxy $minz\"\n",
        vmmc.simbox().min_x(),
        vmmc.simbox().max_x(),
        vmmc.simbox().min_y(),
        vmmc.simbox().max_y()
    )
    .unwrap();
}

fn draw_white_background(vmmc: &Vmmc, dt: &mut DrawTarget, scale: Num) {
    let x_dim = (vmmc.simbox().max_x() * 2.0 * scale) as f32;
    let y_dim = (vmmc.simbox().max_y() * 2.0 * scale) as f32;

    // get white color
    let source = rgba_solid(0xff, 0xff, 0xff, 0xff);

    let mut pb = PathBuilder::new();
    pb.move_to(0.0, 0.0);
    pb.line_to(x_dim, 0.0);
    pb.line_to(x_dim, y_dim);
    pb.line_to(0.0, y_dim);
    pb.line_to(0.0, 0.0);
    let draw_path = pb.finish();
    dt.fill(&draw_path, &source, &DrawOptions::new());
}

fn rgba_solid(r: u8, g: u8, b: u8, a: u8) -> Source<'static> {
    Source::Solid(SolidSource { r, g, b, a })
}

// fn color_polygon(vmmc: &Vmmc, polygon: &Polygon, dt: &mut DrawTarget, scale: Num) {
//     let x_off = vmmc.simbox().max_x() * scale;
//     let y_off = vmmc.simbox().max_y() * scale;

//     let source = match polygon.vertices().len() {
//         0 | 1 | 2 => panic!("Detected \"polygon\" with less than 3 sides!"),
//         3 => rgba_solid(0xff, 0xff, 0, 0xff), // yellow
//         4 => rgba_solid(0, 0, 0xff, 0xff),    // blue
//         5 => rgba_solid(0xff, 0, 0, 0xff),    // red
//         6 => rgba_solid(0x0, 0xff, 0, 0xff),  // green
//         _ => rgba_solid(0xff, 0, 0xff, 0xff), // purple
//     };

//     let mut pb = PathBuilder::new();
//     let mut pos = vmmc.particle(polygon.vertices()[0]).pos();
//     pb.move_to(
//         (pos.x() * scale + x_off) as f32,
//         (pos.y() * scale + y_off) as f32,
//     );

//     for (src_p, dst_p) in polygon.edge_iter() {
//         let pos1 = vmmc.particle(dst_p).pos();
//         let sep = vmmc.simbox().sep_in_box(pos1, pos); // pos -> new_pos
//         let scaled_sep = sep.scalar_mul(scale);
//         let final_pos = pos + sep.scalar_mul(scale);
//         pb.line_to(final_pos.x() as f32, final_pos.y() as f32);
//         pos = final_pos;
//     }

//     let draw_path = pb.finish();
//     dt.fill(&draw_path, &source, &DrawOptions::new());
// }

// fn color_polygons(vmmc: &Vmmc, polygons: &[Polygon], dt: &mut DrawTarget, scale: Num) {
//     for polygon in polygons.iter() {
//         color_polygon(vmmc, polygon, dt, scale);
//     }
// }

fn render_particles(vmmc: &Vmmc, dt: &mut DrawTarget, scale: Num) {
    let x_off = vmmc.simbox().max_x() * scale;
    let y_off = vmmc.simbox().max_y() * scale;
    let source = rgba_solid(0, 0, 0, 0xff);

    let style = StrokeStyle::default();
    let draw_options = DrawOptions::new();

    for p in vmmc.particles().iter() {
        let mut pb = PathBuilder::new();
        let x = (p.pos().x() * scale + x_off) as f32;
        let y = (p.pos().y() * scale + y_off) as f32;
        // pb.move_to(x,y);
        pb.arc(x, y, (scale / 2.0) as f32, 0.0, 2.0 * PI);
        pb.close();
        let draw_path = pb.finish();
        dt.stroke(&draw_path, &source, &style, &draw_options);

        // draw patches
        for (patch_idx, patch) in vmmc.simbox().morphology(p).patches().iter().enumerate() {
            let mut pb = PathBuilder::new();
            let patch_center = vmmc.simbox().patch_center(p, patch_idx);
            let patch_x = (patch_center.x() * scale + x_off) as f32;
            let patch_y = (patch_center.y() * scale + y_off) as f32;
            pb.arc(
                patch_x,
                patch_y,
                (scale * patch.radius()) as f32,
                0.0,
                2.0 * PI,
            );
            pb.close();
            let draw_path = pb.finish();
            dt.stroke(&draw_path, &source, &style, &draw_options);
        }
    }
}

fn render_interactions(vmmc: &Vmmc, dt: &mut DrawTarget, scale: Num) {
    let x_off = vmmc.simbox().max_x() * scale;
    let y_off = vmmc.simbox().max_y() * scale;
    let source = rgba_solid(0, 0, 0, 0xff);

    let style = StrokeStyle::default();
    let draw_options = DrawOptions::new();

    for p0 in vmmc.particles().iter() {
        let interactions = vmmc.potential().determine_interactions(vmmc.simbox(), p0);
        for &neighbor_id in interactions.iter() {
            let p1 = vmmc.particle(neighbor_id);
            // vector from p0 -> p1
            let sep = vmmc.simbox().sep_in_box(p1.pos(), p0.pos());

            let mut pb = PathBuilder::new();
            pb.move_to(
                (p0.pos().x() * scale + x_off) as f32,
                (p0.pos().y() * scale + y_off) as f32,
            );
            let pos1 = p0.pos() + sep;
            pb.line_to(
                (pos1.x() * scale + x_off) as f32,
                (pos1.y() * scale + y_off) as f32,
            );
            pb.close();
            let draw_path = pb.finish();

            dt.stroke(&draw_path, &source, &style, &draw_options);
        }
    }
}

pub fn write_geometry_png(vmmc: &Vmmc, pathname: &str) {
    // need to adjust from 0-centered to all positive coordinates
    let scale = 100.0;
    let x_dim = vmmc.simbox().max_x() * 2.0 * scale;
    let y_dim = vmmc.simbox().max_y() * 2.0 * scale;
    let mut dt = DrawTarget::new(x_dim.ceil() as i32, y_dim.ceil() as i32);

    draw_white_background(vmmc, &mut dt, scale);
    render_particles(vmmc, &mut dt, scale);
    render_interactions(vmmc, &mut dt, scale);

    dt.write_png(pathname).unwrap();
}

fn try_delete(p: String) {
    if std::path::Path::new(&p).is_file() {
        fs::remove_file(p).expect("unreachable")
    }
}

pub fn clear_out_files(config: &VmmcConfig) -> anyhow::Result<()> {
    try_delete(config.toml());
    try_delete(config.vmd());
    try_delete(config.geometry());
    try_delete(config.trajectory());
    Ok(())
}
