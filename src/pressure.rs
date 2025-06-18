use crate::{particle::{Particles,IsParticle,Particle,ParticleId}, position::{DimVec, Position}, simbox::SimBox};
use crate::consts::{MAX_PARTICLES_PER_CELL, PARTICLE_DIAMETER, PARTICLE_RADIUS};

pub fn actual_overlap(p0: Position,p1: Position,simbox: &SimBox) -> bool {
    let diff = simbox.sep_in_box(p0, p1);
    let dist = diff.l2_norm_sqd().sqrt();

    if dist < PARTICLE_DIAMETER {
        return true;
    } 
    false
}

pub fn PressureAtScalingFactor(sim: &SimBox, scalingFactor: f64 ) -> f64 {
    let mut p_vec: Vec<Particle> = vec![];
    
    for old_particle in sim.particles().iter() {
        let new_pos = Position::new([old_particle.pos().x()*(1.0-scalingFactor),old_particle.pos().y()*(1.0-scalingFactor)]);
        p_vec.push(Particle::new(old_particle.id(),new_pos,old_particle.or(),old_particle.shape_id()));
        
    }

    let new_dimensions = DimVec::new([2.0*sim.max_x()*(1.0-scalingFactor),2.0*sim.max_y()*(1.0-scalingFactor)]);
    let new_cell_dimensions = DimVec::new([sim.cell_dimensions().x()*(1.0-scalingFactor),sim.cell_dimensions().y()*(1.0-scalingFactor)]);
    let new_sim_box = SimBox::new(new_dimensions,sim.cells_per_axis(),new_cell_dimensions,Particles::new(p_vec),sim.shapes().to_vec());

    //count overllaps in the new sim box
    let mut N = 0;
    for p in new_sim_box.particles().iter() {
        for p2 in new_sim_box.particles().iter(){
            if p.id() != p2.id() {
                if actual_overlap(p.pos(),p2.pos(),sim) {
                    N += 1;
                }
            }
        }
        //if new_sim_box.overlaps(p) {
        //    N += 1;
        //}
    }
    //println!("{}",new_sim_box.particles().iter().any(|p| new_sim_box.overlaps(p)));
    //println!("overlaps: {}",N);
    return (N as f64)/(2.0*scalingFactor*2.0*(sim.max_x() as f64)*2.0*(sim.max_y() as f64));

}

pub fn Pressure(sim: &SimBox) -> f64 {
    let maxScalingFactor = 0.1;
    let iterations = 10;
    let mut totalPressure = 0.0;
    for i in 0..iterations {
        let scalingFactor = ((i as f64)+1.0)*(maxScalingFactor/(iterations as f64));
        totalPressure += PressureAtScalingFactor(sim, scalingFactor) ;
    } 
    return totalPressure/(iterations as f64);
}

pub fn SimboxAtVolume(sim: &SimBox, volume: DimVec ) -> SimBox {
    let xScalingFactor = volume.x()/(sim.max_x()*2.0);
    let yScalingFactor = volume.y()/(sim.max_y()*2.0);
    
    let mut p_vec: Vec<Particle> = vec![];
    for old_particle in sim.particles().iter() {
        let new_pos = Position::new([old_particle.pos().x()*(xScalingFactor),old_particle.pos().y()*(yScalingFactor)]);
        p_vec.push(Particle::new(old_particle.id(),new_pos,old_particle.or(),old_particle.shape_id()));
    }

    let new_cell_dimensions = DimVec::new([sim.cell_dimensions().x()*(xScalingFactor),sim.cell_dimensions().y()*(1.0-yScalingFactor)]);
    let new_sim_box = SimBox::new(volume,sim.cells_per_axis(),new_cell_dimensions,Particles::new(p_vec),sim.shapes().to_vec());
    return new_sim_box;
}


pub fn Volume(pressure: f64, sim: &SimBox ) -> DimVec {
    //println!("volume started again");
    let mut currentPressure = Pressure(sim);
    let mut dPdX = 1.0;//some random initial change to start
    let mut changeRatio = 1.0;
    //let mut sign = 1.0;
    let discrepancy = 0.01;
    if currentPressure < pressure {
        changeRatio = -1.0;
    }
    let mut accept = false;
    let mut currentVolume = DimVec::new([2.0*sim.max_x()*(1.0),2.0*sim.max_y()*(1.0)]);
    while !accept {
        if currentPressure <= pressure*(1.0 + discrepancy) && currentPressure >= pressure*(1.0 - discrepancy)  {
            accept = true; 
        } else {
            let change = dPdX*changeRatio;
            //currentVolume.set_x(currentVolume.x() + dPdX);
            //currentVolume.y = currentVolume.y + dPdX;
            
            //currentVolume.translate_by(dPdX,dPdX);
            //^ doesnt work at all.
            currentVolume = DimVec::new([currentVolume.x()+change,currentVolume.y()+change]);

            //println!("current volume: {},{}",currentVolume.x(),currentVolume.y());
            let newSim = SimboxAtVolume(sim,currentVolume);
            let newPressure = Pressure(&newSim);

            changeRatio = (newPressure - pressure)/pressure.max(newPressure);
            //gradient descent: x = x - dP/dx
            //dPdX = (newPressure - currentPressure)/dPdX; this does not make sense
            //dPdX = (pressure - newPressure)/dPdX;
            
            
            
            //println!("last pressure {} current pressure {}, dPdX {}, change ratio {}",currentPressure,newPressure,dPdX,changeRatio);
            currentPressure = newPressure;
            //println!("last pressure {} current pressure {}, dPdX {}",currentPressure,dPdX);
        }


    }

    return currentVolume;
}

pub fn simboxAtPressure(pressure: f64, sim: &SimBox) -> SimBox {
    let newVolume = Volume(pressure,sim);
    return SimboxAtVolume(sim,newVolume);
}
