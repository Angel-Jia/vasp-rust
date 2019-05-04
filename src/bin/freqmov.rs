extern crate vasp_rust;

use std::env;
use std::fs::File;
use std::io::BufReader;
use std::io::prelude::*;

fn read_freq(file_name: &str) -> (Vec<(f64, f64, f64)>, Vec<(f64, f64, f64)>){
    let f = File::open(file_name).expect(&format!("can not open file: {}", file_name));
    let input_file = BufReader::new(f).lines().map(|s| s.unwrap()).collect::<Vec<String>>();
    
    let mut coord: Vec<(f64, f64, f64)> = Vec::new();
    let mut diff: Vec<(f64, f64, f64)> = Vec::new();

    for line in input_file[2..].iter(){
        let line = line.trim();
        if line.len() == 0 {break;}
        let tmp = line.split_whitespace().map(|num| num.parse::<f64>().unwrap()).collect::<Vec<f64>>();
        coord.push((tmp[0], tmp[1], tmp[2]));
        diff.push((tmp[3], tmp[4], tmp[5]));
    }

    (coord, diff)
}

fn coord_mul(vector1: &Vec<(f64, f64, f64)>, vector2: &Vec<(f64, f64, f64)>, scale: f64, ret: &mut Vec<(f64, f64, f64)>){
    for (v1, v2) in vector1.iter().zip(vector2.iter()){
        ret.push((v1.0 + v2.0 * scale, v1.1 + v2.1 * scale, v1.2 + v2.2 * scale));
    }
}

fn main(){
    let mut args: Vec<String> = env::args().collect();
    if args.len() < 4{
        println!("\nUsage: {} freq_file1 freq_file2 .. frame_num scale", args[0]);
        println!("Please try again!\n");
    }
    let scale_s = args.pop().unwrap();
    let scale: f64 = scale_s.parse::<f64>().expect(&format!("can not recgnize ratio: {}", scale_s));

    let frame_num_s = args.pop().unwrap();
    let frame_num: u32 = frame_num_s.parse::<u32>().expect(&format!("can not recgnize ratio: {}", frame_num_s));


    println!("\n################ This script makes animation of vibration ################\n");

    let mut poscar = vasp_rust::Data::new();
    poscar.read_vasp("POSCAR");

    let total_atoms: u64 = poscar.num_atoms.iter().sum();

    for freq_name in args[1..].iter(){
        println!("                            Processing {}", freq_name);
        let (coord, diff) = read_freq(freq_name);

        if coord.len() != total_atoms as usize {panic!(format!("atoms number({}) in freq file differ from that({}) in POSCAR",
            coord.len(), total_atoms))}
        
        let mut ret: Vec<(f64, f64, f64)> = Vec::with_capacity(total_atoms as usize * (frame_num as usize * 4 + 4));
        
        //0
        ret.extend(coord.clone());

        //(0, 1)
        let scale = 1.0 * scale / (frame_num as f64);
        for s in 1..frame_num{
            coord_mul(&coord, &diff, s as f64 * scale, &mut ret);
        }

        //1
        coord_mul(&coord, &diff, frame_num as f64 * scale, &mut ret);

        //(1 → -1]
        for s in (-(frame_num as i64)..frame_num as i64).rev(){
            coord_mul(&coord, &diff, s as f64 * scale, &mut ret);
        }

        //(-1, 0)
        for s in -((frame_num - 1) as i64)..0{
            coord_mul(&coord, &diff, s as f64 * scale, &mut ret);
        }

        poscar.coordinates = ret;
        poscar.write_xyz(&format!("{}.xyz", freq_name));
    }

    println!("\n                ------------------ Done ------------------\n")
}
