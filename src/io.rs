use std::{
    fs::File,
    io::{BufRead, BufReader, Write},
};

use crate::{particle::IsParticle, position::DimVec, vmmc::Vmmc};

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
