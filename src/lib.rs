mod utils;

use wasm_bindgen::prelude::*;
use js_sys::{Float64Array, Uint32Array};

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

#[derive(Clone, Copy)]
pub enum IntegrationMethod {
    Euler,
    RungeKutta4,
}

#[derive(Clone, Copy)]
pub enum AttractorType {
    Lorenz,
    Rossler,
    Thomas,
    Chua,
}

#[wasm_bindgen]
pub struct ChaoticAttractor {
    x: f64,
    y: f64,
    z: f64,
    // Parameters - their meaning depends on attractor type
    param_a: f64,
    param_b: f64,
    param_c: f64,
    dt: f64,
    points: Vec<Point3D>,
    max_points: usize,
    integration_method: IntegrationMethod,
    buffer_start: usize,
    buffer_count: usize,
    attractor_type: AttractorType,
    // For chaos analysis
    lyapunov_sum: f64,
    lyapunov_steps: u32,
    poincare_plane_z: f64,
    poincare_points: Vec<Point3D>,
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
impl ChaoticAttractor {
    #[wasm_bindgen(constructor)]
    pub fn new() -> ChaoticAttractor {
        utils::set_panic_hook();
        
        let max_points = 10000;
        ChaoticAttractor {
            x: 0.1,
            y: 0.0,
            z: 0.0,
            param_a: 10.0,  // sigma for Lorenz
            param_b: 28.0,  // rho for Lorenz
            param_c: 8.0 / 3.0,  // beta for Lorenz
            dt: 0.01,
            points: Vec::with_capacity(max_points),
            max_points,
            integration_method: IntegrationMethod::Euler,
            buffer_start: 0,
            buffer_count: 0,
            attractor_type: AttractorType::Lorenz,
            lyapunov_sum: 0.0,
            lyapunov_steps: 0,
            poincare_plane_z: 27.0,
            poincare_points: Vec::new(),
        }
    }
    
    #[wasm_bindgen]
    pub fn set_parameters(&mut self, a: f64, b: f64, c: f64) {
        self.param_a = a;
        self.param_b = b;
        self.param_c = c;
    }
    
    #[wasm_bindgen]
    pub fn set_attractor_type(&mut self, attractor_type: u32) {
        self.attractor_type = match attractor_type {
            1 => AttractorType::Rossler,
            2 => AttractorType::Thomas,
            3 => AttractorType::Chua,
            _ => AttractorType::Lorenz,
        };
        
        // Set default parameters for each attractor
        match self.attractor_type {
            AttractorType::Lorenz => {
                self.param_a = 10.0;   // sigma
                self.param_b = 28.0;   // rho  
                self.param_c = 8.0/3.0; // beta
                self.x = 0.1;
                self.y = 0.0;
                self.z = 0.0;
            },
            AttractorType::Rossler => {
                self.param_a = 0.2;    // a
                self.param_b = 0.2;    // b
                self.param_c = 5.7;    // c
                self.x = 1.0;
                self.y = 1.0;
                self.z = 1.0;
            },
            AttractorType::Thomas => {
                self.param_a = 0.208186; // b (Thomas parameter)
                self.param_b = 0.0;      // unused
                self.param_c = 0.0;      // unused
                self.x = 0.1;
                self.y = 0.0;
                self.z = 0.0;
            },
            AttractorType::Chua => {
                self.param_a = 15.6;   // alpha
                self.param_b = 28.0;   // beta
                self.param_c = -1.143; // m0
                self.x = 0.7;
                self.y = 0.0;
                self.z = 0.0;
            },
        }
    }
    
    #[wasm_bindgen]
    pub fn set_initial_conditions(&mut self, x: f64, y: f64, z: f64) {
        self.x = x;
        self.y = y;
        self.z = z;
        self.buffer_start = 0;
        self.buffer_count = 0;
    }
    
    #[wasm_bindgen]
    pub fn set_time_step(&mut self, dt: f64) {
        self.dt = dt;
    }
    
    #[wasm_bindgen]
    pub fn set_max_points(&mut self, max_points: usize) {
        self.max_points = max_points;
        if self.buffer_count > max_points {
            let excess = self.buffer_count - max_points;
            self.buffer_start = (self.buffer_start + excess) % self.points.len();
            self.buffer_count = max_points;
        }
        
        if self.points.capacity() < max_points {
            self.points.reserve(max_points - self.points.capacity());
        }
    }
    
    #[wasm_bindgen]
    pub fn set_integration_method(&mut self, method: u32) {
        self.integration_method = match method {
            1 => IntegrationMethod::RungeKutta4,
            _ => IntegrationMethod::Euler,
        };
    }
    
    fn compute_derivatives(&self, x: f64, y: f64, z: f64) -> (f64, f64, f64) {
        match self.attractor_type {
            AttractorType::Lorenz => {
                let dx = self.param_a * (y - x);
                let dy = x * (self.param_b - z) - y;
                let dz = x * y - self.param_c * z;
                (dx, dy, dz)
            },
            AttractorType::Rossler => {
                let dx = -y - z;
                let dy = x + self.param_a * y;
                let dz = self.param_b + z * (x - self.param_c);
                (dx, dy, dz)
            },
            AttractorType::Thomas => {
                let dx = -self.param_a * x + y.sin();
                let dy = -self.param_a * y + z.sin();  
                let dz = -self.param_a * z + x.sin();
                (dx, dy, dz)
            },
            AttractorType::Chua => {
                // Chua's circuit with piecewise linear nonlinearity
                let m1 = -1.143;  // m1
                let m0 = self.param_c; // m0  
                let f_x = if x.abs() <= 1.0 {
                    m1 * x
                } else {
                    m0 * x + (m1 - m0) * x.signum()
                };
                
                let dx = self.param_a * (y - x - f_x);
                let dy = x - y + z;
                let dz = -self.param_b * y;
                (dx, dy, dz)
            }
        }
    }
    
    fn compute_jacobian(&self, x: f64, y: f64, z: f64) -> [[f64; 3]; 3] {
        match self.attractor_type {
            AttractorType::Lorenz => [
                [-self.param_a, self.param_a, 0.0],
                [self.param_b - z, -1.0, -x],
                [y, x, -self.param_c]
            ],
            AttractorType::Rossler => [
                [0.0, -1.0, -1.0],
                [1.0, self.param_a, 0.0],
                [z, 0.0, x - self.param_c]
            ],
            AttractorType::Thomas => [
                [-self.param_a, y.cos(), 0.0],
                [0.0, -self.param_a, z.cos()],
                [x.cos(), 0.0, -self.param_a]
            ],
            AttractorType::Chua => {
                let df_dx = if x.abs() <= 1.0 {
                    -1.143
                } else {
                    self.param_c
                };
                [
                    [self.param_a * (-1.0 - df_dx), self.param_a, 0.0],
                    [1.0, -1.0, 1.0],
                    [0.0, -self.param_b, 0.0]
                ]
            }
        }
    }
    
    fn matrix_norm_3x3(m: [[f64; 3]; 3]) -> f64 {
        let mut sum = 0.0;
        for i in 0..3 {
            for j in 0..3 {
                sum += m[i][j] * m[i][j];
            }
        }
        sum.sqrt()
    }
    
    #[wasm_bindgen]
    pub fn step(&mut self) {
        match self.integration_method {
            IntegrationMethod::Euler => {
                let (dx, dy, dz) = self.compute_derivatives(self.x, self.y, self.z);
                self.x += dx * self.dt;
                self.y += dy * self.dt;
                self.z += dz * self.dt;
            },
            IntegrationMethod::RungeKutta4 => {
                let h = self.dt;
                let (k1x, k1y, k1z) = self.compute_derivatives(self.x, self.y, self.z);
                
                let (k2x, k2y, k2z) = self.compute_derivatives(
                    self.x + 0.5 * h * k1x,
                    self.y + 0.5 * h * k1y,
                    self.z + 0.5 * h * k1z
                );
                
                let (k3x, k3y, k3z) = self.compute_derivatives(
                    self.x + 0.5 * h * k2x,
                    self.y + 0.5 * h * k2y,
                    self.z + 0.5 * h * k2z
                );
                
                let (k4x, k4y, k4z) = self.compute_derivatives(
                    self.x + h * k3x,
                    self.y + h * k3y,
                    self.z + h * k3z
                );
                
                self.x += h * (k1x + 2.0 * k2x + 2.0 * k3x + k4x) / 6.0;
                self.y += h * (k1y + 2.0 * k2y + 2.0 * k3y + k4y) / 6.0;
                self.z += h * (k1z + 2.0 * k2z + 2.0 * k3z + k4z) / 6.0;
            }
        }
        
        // Store previous z for Poincaré analysis
        let prev_z = if self.buffer_count > 0 {
            self.get_point(self.buffer_count - 1).map_or(self.z, |p| p.z)
        } else {
            self.z
        };
        
        let new_point = Point3D {
            x: self.x,
            y: self.y,
            z: self.z,
        };
        
        if self.points.len() < self.max_points {
            self.points.push(new_point);
            self.buffer_count += 1;
        } else {
            let insert_index = (self.buffer_start + self.buffer_count) % self.points.len();
            if insert_index < self.points.len() {
                self.points[insert_index] = new_point;
            } else {
                self.points.push(new_point);
            }
            
            if self.buffer_count < self.max_points {
                self.buffer_count += 1;
            } else {
                self.buffer_start = (self.buffer_start + 1) % self.points.len();
            }
        }
        
        // Perform chaos analysis
        self.calculate_lyapunov_step();
        self.check_poincare_intersection(prev_z, self.z);
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
        self.buffer_count
    }
    
    #[wasm_bindgen]
    pub fn get_points_bulk(&self) -> Float64Array {
        let mut data = Vec::with_capacity(self.buffer_count * 3);
        
        for i in 0..self.buffer_count {
            let idx = (self.buffer_start + i) % self.points.len();
            if idx < self.points.len() {
                let point = &self.points[idx];
                data.push(point.x);
                data.push(point.y);
                data.push(point.z);
            }
        }
        
        Float64Array::from(&data[..])
    }
    
    #[wasm_bindgen]
    pub fn get_memory_stats(&self) -> Uint32Array {
        let stats = vec![
            self.points.len() as u32,
            self.points.capacity() as u32, 
            self.buffer_count as u32,
            self.max_points as u32,
        ];
        Uint32Array::from(&stats[..])
    }
    
    #[wasm_bindgen]
    pub fn get_point(&self, index: usize) -> Option<Point3D> {
        if index >= self.buffer_count {
            return None;
        }
        
        let actual_index = (self.buffer_start + index) % self.points.len();
        self.points.get(actual_index).map(|p| Point3D {
            x: p.x,
            y: p.y,
            z: p.z,
        })
    }
    
    #[wasm_bindgen]
    pub fn clear_points(&mut self) {
        self.buffer_start = 0;
        self.buffer_count = 0;
    }
    
    #[wasm_bindgen]
    pub fn get_x_range(&self) -> Vec<f64> {
        if self.buffer_count == 0 {
            return vec![-20.0, 20.0];
        }
        
        let first_idx = self.buffer_start % self.points.len();
        let mut min_x = self.points[first_idx].x;
        let mut max_x = self.points[first_idx].x;
        
        for i in 0..self.buffer_count {
            let idx = (self.buffer_start + i) % self.points.len();
            if idx < self.points.len() {
                let point = &self.points[idx];
                if point.x < min_x { min_x = point.x; }
                if point.x > max_x { max_x = point.x; }
            }
        }
        
        vec![min_x, max_x]
    }
    
    #[wasm_bindgen]
    pub fn get_y_range(&self) -> Vec<f64> {
        if self.buffer_count == 0 {
            return vec![-30.0, 30.0];
        }
        
        let first_idx = self.buffer_start % self.points.len();
        let mut min_y = self.points[first_idx].y;
        let mut max_y = self.points[first_idx].y;
        
        for i in 0..self.buffer_count {
            let idx = (self.buffer_start + i) % self.points.len();
            if idx < self.points.len() {
                let point = &self.points[idx];
                if point.y < min_y { min_y = point.y; }
                if point.y > max_y { max_y = point.y; }
            }
        }
        
        vec![min_y, max_y]
    }
    
    #[wasm_bindgen]
    pub fn get_z_range(&self) -> Vec<f64> {
        if self.buffer_count == 0 {
            return vec![0.0, 50.0];
        }
        
        let first_idx = self.buffer_start % self.points.len();
        let mut min_z = self.points[first_idx].z;
        let mut max_z = self.points[first_idx].z;
        
        for i in 0..self.buffer_count {
            let idx = (self.buffer_start + i) % self.points.len();
            if idx < self.points.len() {
                let point = &self.points[idx];
                if point.z < min_z { min_z = point.z; }
                if point.z > max_z { max_z = point.z; }
            }
        }
        
        vec![min_z, max_z]
    }
    
    #[wasm_bindgen]
    pub fn get_attractor_type(&self) -> u32 {
        match self.attractor_type {
            AttractorType::Lorenz => 0,
            AttractorType::Rossler => 1,
            AttractorType::Thomas => 2,
            AttractorType::Chua => 3,
        }
    }
    
    #[wasm_bindgen]
    pub fn get_parameters(&self) -> Vec<f64> {
        vec![self.param_a, self.param_b, self.param_c]
    }
    
    #[wasm_bindgen]
    pub fn get_attractor_info(&self) -> String {
        match self.attractor_type {
            AttractorType::Lorenz => "Lorenz|σ (sigma)|ρ (rho)|β (beta)|0.1,20|0.1,50|0.1,5".to_string(),
            AttractorType::Rossler => "Rössler|a|b|c|0.01,1|0.01,1|0.1,20".to_string(),
            AttractorType::Thomas => "Thomas|b|unused|unused|0.01,1|0,0|0,0".to_string(),
            AttractorType::Chua => "Chua|α (alpha)|β (beta)|m₀|1,50|1,50|-5,5".to_string(),
        }
    }
    
    #[wasm_bindgen]
    pub fn calculate_lyapunov_step(&mut self) -> f64 {
        let jacobian = self.compute_jacobian(self.x, self.y, self.z);
        
        // Use Frobenius norm of Jacobian as approximation for largest eigenvalue
        let norm = Self::matrix_norm_3x3(jacobian);
        let lyap_contrib = norm.ln();
        
        self.lyapunov_sum += lyap_contrib;
        self.lyapunov_steps += 1;
        
        lyap_contrib
    }
    
    #[wasm_bindgen]
    pub fn get_lyapunov_exponent(&self) -> f64 {
        if self.lyapunov_steps == 0 {
            return 0.0;
        }
        self.lyapunov_sum / (self.lyapunov_steps as f64 * self.dt)
    }
    
    #[wasm_bindgen]
    pub fn reset_lyapunov(&mut self) {
        self.lyapunov_sum = 0.0;
        self.lyapunov_steps = 0;
    }
    
    #[wasm_bindgen]
    pub fn set_poincare_plane(&mut self, z_value: f64) {
        self.poincare_plane_z = z_value;
        self.poincare_points.clear();
    }
    
    fn check_poincare_intersection(&mut self, prev_z: f64, current_z: f64) {
        let z_plane = self.poincare_plane_z;
        
        // Check if trajectory crossed the plane (upward crossing only)
        if prev_z < z_plane && current_z >= z_plane {
            // Linear interpolation to find intersection point
            let t = (z_plane - prev_z) / (current_z - prev_z);
            
            // Get the previous point for interpolation
            if self.buffer_count >= 2 {
                if let Some(prev_point) = self.get_point(self.buffer_count - 2) {
                    let intersect_x = prev_point.x + t * (self.x - prev_point.x);
                    let intersect_y = prev_point.y + t * (self.y - prev_point.y);
                    
                    self.poincare_points.push(Point3D {
                        x: intersect_x,
                        y: intersect_y,
                        z: z_plane,
                    });
                    
                    // Keep only recent Poincaré points
                    if self.poincare_points.len() > 1000 {
                        self.poincare_points.remove(0);
                    }
                }
            }
        }
    }
    
    #[wasm_bindgen]
    pub fn get_poincare_points(&self) -> Float64Array {
        let mut data = Vec::with_capacity(self.poincare_points.len() * 2);
        for point in &self.poincare_points {
            data.push(point.x);
            data.push(point.y);
        }
        Float64Array::from(&data[..])
    }
    
    #[wasm_bindgen]
    pub fn get_poincare_count(&self) -> usize {
        self.poincare_points.len()
    }
    
    #[wasm_bindgen]
    pub fn clear_poincare_points(&mut self) {
        self.poincare_points.clear();
    }
}