use std::{
    fs::File,
    io::{BufRead, BufReader, Write},
};

// use raqote::{PathBuilder, DrawTarget, Source, SolidSource, StrokeStyle, DrawOptions};
use raqote::*;

use crate::{particle::IsParticle, position::DimVec, simbox, vmmc::Vmmc};

pub struct XYZWriter {
    file: File,
}

pub struct XYZReader {
    rdr: BufReader<File>,
}

impl XYZReader {
    pub fn new(p: &str) -> Self {
        let file = File::open(p).unwrap();
        let rdr = BufReader::new(file);
        Self { rdr }
    }
}

pub fn read_xyz_snapshot(path: &str) -> (Vec<DimVec>, Vec<DimVec>) {
    let mut positions = Vec::new();
    let mut orientations = Vec::new();

    let file = File::open(path).unwrap();
    let rdr = BufReader::new(file);

    for line in rdr.lines() {
        let line = line.unwrap();
        let parts: Vec<&str> = line.split_whitespace().collect();
        let pos_x = parts[0].parse::<f64>().unwrap();
        let pos_y = parts[1].parse::<f64>().unwrap();
        let or_x = parts[2].parse::<f64>().unwrap();
        let or_y = parts[3].parse::<f64>().unwrap();
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
        writeln!(self.file, "{:?}\n", vmmc.particles().len()).unwrap();
        for p in vmmc.particles() {
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

fn draw_white_background(vmmc: &Vmmc, dt: &mut DrawTarget, scale: f64) {
    let x_dim = (vmmc.simbox().max_x() * 2.0 * scale) as f32;
    let y_dim = (vmmc.simbox().max_y() * 2.0 * scale) as f32;

    let source = Source::Solid(SolidSource {
        r: 0xff,
        g: 0xff,
        b: 0xff,
        a: 0xff,
    });
    // let style = StrokeStyle::default();
    // let draw_options = DrawOptions::new();

    let mut pb = PathBuilder::new();
    pb.move_to(0.0, 0.0);
    pb.line_to(x_dim, 0.0);
    pb.line_to(x_dim, y_dim);
    pb.line_to(0.0, y_dim);
    pb.line_to(0.0, 0.0);
    let draw_path = pb.finish();
    // dt.stroke(&draw_path, &source, &style, &draw_options);
    dt.fill(&draw_path, &source, &DrawOptions::new());
}

pub fn write_geometry_png(vmmc: &Vmmc, pathname: &str) {
    // need to adjust from 0-centered to all positive coordinates
    let scale = 20.0;
    let x_dim = vmmc.simbox().max_x() * 2.0 * scale;
    let y_dim = vmmc.simbox().max_y() * 2.0 * scale;
    let mut dt = DrawTarget::new(x_dim.ceil() as i32, y_dim.ceil() as i32);
    let x_off = vmmc.simbox().max_x() * scale;
    let y_off = vmmc.simbox().max_y() * scale;

    draw_white_background(vmmc, &mut dt, scale);

    let source = Source::Solid(SolidSource {
        r: 0x0,
        g: 0x0,
        b: 0x0,
        a: 0xff,
    });
    let style = StrokeStyle {
        cap: LineCap::Round,
        join: LineJoin::Round,
        width: 1.,
        miter_limit: 10.,
        dash_array: Vec::new(),
        dash_offset: 0.,
    };

    let style = StrokeStyle::default();
    let draw_options = DrawOptions::new();

    for p0 in vmmc.particles() {
        let interactions = vmmc.potential().determine_interactions(vmmc.simbox(), p0);
        for &neighbor_id in interactions.iter() {
            let p1 = vmmc.particle(neighbor_id);

            let mut pb = PathBuilder::new();
            pb.move_to(
                (p0.pos().x() * scale + x_off) as f32,
                (p0.pos().y() * scale + y_off) as f32,
            );
            pb.line_to(
                (p1.pos().x() * scale + x_off) as f32,
                (p1.pos().y() * scale + y_off) as f32,
            );
            pb.close();
            let draw_path = pb.finish();

            dt.stroke(&draw_path, &source, &style, &draw_options);
        }
    }
    dt.write_png(pathname).unwrap();
}

pub fn draw_example() {
    let mut dt = DrawTarget::new(400, 400);

    let mut pb = PathBuilder::new();
    pb.move_to(100., 10.);
    pb.cubic_to(150., 40., 175., 0., 200., 10.);
    pb.quad_to(120., 100., 80., 200.);
    pb.quad_to(150., 180., 300., 300.);
    pb.close();
    let path = pb.finish();

    let gradient = Source::new_radial_gradient(
        Gradient {
            stops: vec![
                GradientStop {
                    position: 0.2,
                    color: Color::new(0xff, 0, 0xff, 0),
                },
                GradientStop {
                    position: 0.8,
                    color: Color::new(0xff, 0xff, 0xff, 0xff),
                },
                GradientStop {
                    position: 1.,
                    color: Color::new(0xff, 0xff, 0, 0xff),
                },
            ],
        },
        Point::new(150., 150.),
        128.,
        Spread::Pad,
    );
    dt.fill(&path, &gradient, &DrawOptions::new());

    let mut pb = PathBuilder::new();
    pb.move_to(100., 100.);
    pb.line_to(300., 300.);
    pb.line_to(200., 300.);
    let path = pb.finish();

    dt.stroke(
        &path,
        &Source::Solid(SolidSource {
            r: 0x0,
            g: 0x0,
            b: 0x80,
            a: 0x80,
        }),
        &StrokeStyle {
            cap: LineCap::Round,
            join: LineJoin::Round,
            width: 10.,
            miter_limit: 2.,
            dash_array: vec![10., 18.],
            dash_offset: 16.,
        },
        &DrawOptions::new(),
    );

    dt.write_png("example.png");
}
