mod utils;

use wasm_bindgen::prelude::*;

#[cfg(feature = "wee_alloc")]
#[global_allocator]
static ALLOC: wee_alloc::WeeAlloc = wee_alloc::WeeAlloc::INIT;

#[wasm_bindgen]
extern "C" {
    fn alert(s: &str);
    
    #[wasm_bindgen(js_namespace = console)]
    fn log(s: &str);
}

macro_rules! console_log {
    ($($t:tt)*) => (log(&format_args!($($t)*).to_string()))
}

#[wasm_bindgen]
pub struct LorenzAttractor {
    x: f64,
    y: f64,
    z: f64,
    sigma: f64,
    rho: f64,
    beta: f64,
    dt: f64,
    points: Vec<Point3D>,
    max_points: usize,
}

#[wasm_bindgen]
pub struct Point3D {
    x: f64,
    y: f64,
    z: f64,
}

#[wasm_bindgen]
impl Point3D {
    #[wasm_bindgen(getter)]
    pub fn x(&self) -> f64 {
        self.x
    }
    
    #[wasm_bindgen(getter)]
    pub fn y(&self) -> f64 {
        self.y
    }
    
    #[wasm_bindgen(getter)]
    pub fn z(&self) -> f64 {
        self.z
    }
}

#[wasm_bindgen]
impl LorenzAttractor {
    #[wasm_bindgen(constructor)]
    pub fn new() -> LorenzAttractor {
        utils::set_panic_hook();
        
        LorenzAttractor {
            x: 0.1,
            y: 0.0,
            z: 0.0,
            sigma: 10.0,
            rho: 28.0,
            beta: 8.0 / 3.0,
            dt: 0.01,
            points: Vec::new(),
            max_points: 10000,
        }
    }
    
    #[wasm_bindgen]
    pub fn set_parameters(&mut self, sigma: f64, rho: f64, beta: f64) {
        self.sigma = sigma;
        self.rho = rho;
        self.beta = beta;
    }
    
    #[wasm_bindgen]
    pub fn set_initial_conditions(&mut self, x: f64, y: f64, z: f64) {
        self.x = x;
        self.y = y;
        self.z = z;
        self.points.clear();
    }
    
    #[wasm_bindgen]
    pub fn set_time_step(&mut self, dt: f64) {
        self.dt = dt;
    }
    
    #[wasm_bindgen]
    pub fn set_max_points(&mut self, max_points: usize) {
        self.max_points = max_points;
        if self.points.len() > max_points {
            self.points.drain(0..self.points.len() - max_points);
        }
    }
    
    #[wasm_bindgen]
    pub fn step(&mut self) {
        let dx = self.sigma * (self.y - self.x);
        let dy = self.x * (self.rho - self.z) - self.y;
        let dz = self.x * self.y - self.beta * self.z;
        
        self.x += dx * self.dt;
        self.y += dy * self.dt;
        self.z += dz * self.dt;
        
        self.points.push(Point3D {
            x: self.x,
            y: self.y,
            z: self.z,
        });
        
        if self.points.len() > self.max_points {
            self.points.remove(0);
        }
    }
    
    #[wasm_bindgen]
    pub fn step_multiple(&mut self, steps: usize) {
        for _ in 0..steps {
            self.step();
        }
    }
    
    #[wasm_bindgen]
    pub fn get_current_position(&self) -> Point3D {
        Point3D {
            x: self.x,
            y: self.y,
            z: self.z,
        }
    }
    
    #[wasm_bindgen]
    pub fn get_points_count(&self) -> usize {
        self.points.len()
    }
    
    #[wasm_bindgen]
    pub fn get_point(&self, index: usize) -> Option<Point3D> {
        self.points.get(index).map(|p| Point3D {
            x: p.x,
            y: p.y,
            z: p.z,
        })
    }
    
    #[wasm_bindgen]
    pub fn clear_points(&mut self) {
        self.points.clear();
    }
    
    #[wasm_bindgen]
    pub fn get_x_range(&self) -> Vec<f64> {
        if self.points.is_empty() {
            return vec![-20.0, 20.0];
        }
        
        let mut min_x = self.points[0].x;
        let mut max_x = self.points[0].x;
        
        for point in &self.points {
            if point.x < min_x { min_x = point.x; }
            if point.x > max_x { max_x = point.x; }
        }
        
        vec![min_x, max_x]
    }
    
    #[wasm_bindgen]
    pub fn get_y_range(&self) -> Vec<f64> {
        if self.points.is_empty() {
            return vec![-30.0, 30.0];
        }
        
        let mut min_y = self.points[0].y;
        let mut max_y = self.points[0].y;
        
        for point in &self.points {
            if point.y < min_y { min_y = point.y; }
            if point.y > max_y { max_y = point.y; }
        }
        
        vec![min_y, max_y]
    }
    
    #[wasm_bindgen]
    pub fn get_z_range(&self) -> Vec<f64> {
        if self.points.is_empty() {
            return vec![0.0, 50.0];
        }
        
        let mut min_z = self.points[0].z;
        let mut max_z = self.points[0].z;
        
        for point in &self.points {
            if point.z < min_z { min_z = point.z; }
            if point.z > max_z { max_z = point.z; }
        }
        
        vec![min_z, max_z]
    }
}