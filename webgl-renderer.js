export class WebGLRenderer {
    constructor(canvas) {
        this.canvas = canvas;
        this.gl = canvas.getContext('webgl') || canvas.getContext('experimental-webgl');
        
        if (!this.gl) {
            throw new Error('WebGL not supported');
        }
        
        this.initialize();
    }
    
    initialize() {
        const gl = this.gl;
        
        // Vertex shader for points
        const vertexShaderSource = `
            attribute vec3 a_position;
            attribute float a_age;
            uniform mat4 u_mvpMatrix;
            uniform float u_pointSize;
            varying float v_age;
            
            void main() {
                gl_Position = u_mvpMatrix * vec4(a_position, 1.0);
                gl_PointSize = u_pointSize;
                v_age = a_age;
            }
        `;
        
        // Fragment shader for points
        const fragmentShaderSource = `
            precision mediump float;
            varying float v_age;
            uniform vec3 u_startColor;
            uniform vec3 u_endColor;
            
            void main() {
                // Create circular points
                vec2 coord = gl_PointCoord - vec2(0.5);
                if (length(coord) > 0.5) discard;
                
                // Color interpolation based on age
                vec3 color = mix(u_startColor, u_endColor, v_age);
                float alpha = 0.3 + v_age * 0.7;
                
                gl_FragColor = vec4(color, alpha);
            }
        `;
        
        // Line shader for trails
        const lineVertexShaderSource = `
            attribute vec3 a_position;
            attribute float a_age;
            uniform mat4 u_mvpMatrix;
            varying float v_age;
            
            void main() {
                gl_Position = u_mvpMatrix * vec4(a_position, 1.0);
                v_age = a_age;
            }
        `;
        
        const lineFragmentShaderSource = `
            precision mediump float;
            varying float v_age;
            uniform vec3 u_startColor;
            uniform vec3 u_endColor;
            
            void main() {
                vec3 color = mix(u_startColor, u_endColor, v_age);
                float alpha = 0.1 + v_age * 0.5;
                gl_FragColor = vec4(color, alpha);
            }
        `;
        
        this.pointProgram = this.createProgram(vertexShaderSource, fragmentShaderSource);
        this.lineProgram = this.createProgram(lineVertexShaderSource, lineFragmentShaderSource);
        
        this.setupUniforms();
        this.setupBuffers();
        
        // Enable blending for transparency
        gl.enable(gl.BLEND);
        gl.blendFunc(gl.SRC_ALPHA, gl.ONE_MINUS_SRC_ALPHA);
        gl.enable(gl.DEPTH_TEST);
        
        // Camera settings
        this.camera = {
            distance: 100,
            rotationX: 0.2,
            rotationY: 0,
            autoRotate: false,
            rotationSpeed: 0.01
        };
        
        this.colorScheme = 0; // 0: blue-cyan, 1: rainbow, 2: fire, 3: plasma
    }
    
    createShader(source, type) {
        const gl = this.gl;
        const shader = gl.createShader(type);
        gl.shaderSource(shader, source);
        gl.compileShader(shader);
        
        if (!gl.getShaderParameter(shader, gl.COMPILE_STATUS)) {
            console.error('Shader compilation error:', gl.getShaderInfoLog(shader));
            gl.deleteShader(shader);
            return null;
        }
        
        return shader;
    }
    
    createProgram(vertexSource, fragmentSource) {
        const gl = this.gl;
        const vertexShader = this.createShader(vertexSource, gl.VERTEX_SHADER);
        const fragmentShader = this.createShader(fragmentSource, gl.FRAGMENT_SHADER);
        
        const program = gl.createProgram();
        gl.attachShader(program, vertexShader);
        gl.attachShader(program, fragmentShader);
        gl.linkProgram(program);
        
        if (!gl.getProgramParameter(program, gl.LINK_STATUS)) {
            console.error('Program linking error:', gl.getProgramInfoLog(program));
            gl.deleteProgram(program);
            return null;
        }
        
        return program;
    }
    
    setupUniforms() {
        const gl = this.gl;
        
        // Point program uniforms
        gl.useProgram(this.pointProgram);
        this.pointUniforms = {
            mvpMatrix: gl.getUniformLocation(this.pointProgram, 'u_mvpMatrix'),
            pointSize: gl.getUniformLocation(this.pointProgram, 'u_pointSize'),
            startColor: gl.getUniformLocation(this.pointProgram, 'u_startColor'),
            endColor: gl.getUniformLocation(this.pointProgram, 'u_endColor')
        };
        
        this.pointAttributes = {
            position: gl.getAttribLocation(this.pointProgram, 'a_position'),
            age: gl.getAttribLocation(this.pointProgram, 'a_age')
        };
        
        // Line program uniforms
        gl.useProgram(this.lineProgram);
        this.lineUniforms = {
            mvpMatrix: gl.getUniformLocation(this.lineProgram, 'u_mvpMatrix'),
            startColor: gl.getUniformLocation(this.lineProgram, 'u_startColor'),
            endColor: gl.getUniformLocation(this.lineProgram, 'u_endColor')
        };
        
        this.lineAttributes = {
            position: gl.getAttribLocation(this.lineProgram, 'a_position'),
            age: gl.getAttribLocation(this.lineProgram, 'a_age')
        };
    }
    
    setupBuffers() {
        const gl = this.gl;
        
        this.positionBuffer = gl.createBuffer();
        this.ageBuffer = gl.createBuffer();
        this.indexBuffer = gl.createBuffer();
    }
    
    setColorScheme(scheme) {
        this.colorScheme = scheme;
    }
    
    getColors() {
        const schemes = [
            { start: [0.0, 0.5, 1.0], end: [0.0, 1.0, 1.0] }, // Blue to cyan
            { start: [1.0, 0.0, 0.5], end: [0.0, 1.0, 0.5] }, // Magenta to green
            { start: [1.0, 0.0, 0.0], end: [1.0, 1.0, 0.0] }, // Red to yellow
            { start: [0.5, 0.0, 1.0], end: [1.0, 0.0, 0.5] }  // Purple to magenta
        ];
        
        return schemes[this.colorScheme] || schemes[0];
    }
    
    createMVPMatrix(view) {
        const aspect = this.canvas.width / this.canvas.height;
        
        // Projection matrix
        const fov = Math.PI / 4;
        const near = 0.1;
        const far = 1000.0;
        const f = 1.0 / Math.tan(fov / 2);
        
        const projection = new Float32Array([
            f / aspect, 0, 0, 0,
            0, f, 0, 0,
            0, 0, (far + near) / (near - far), -1,
            0, 0, (2 * far * near) / (near - far), 0
        ]);
        
        // View matrix based on current view mode
        let modelView;
        
        if (view === '3d') {
            const rotX = this.camera.rotationX;
            const rotY = this.camera.rotationY;
            const dist = this.camera.distance;
            
            const cosX = Math.cos(rotX);
            const sinX = Math.sin(rotX);
            const cosY = Math.cos(rotY);
            const sinY = Math.sin(rotY);
            
            modelView = new Float32Array([
                cosY, sinX * sinY, -cosX * sinY, 0,
                0, cosX, sinX, 0,
                sinY, -sinX * cosY, cosX * cosY, 0,
                0, 0, -dist, 1
            ]);
        } else {
            // 2D projections
            const scale = 0.1;
            modelView = new Float32Array([
                scale, 0, 0, 0,
                0, scale, 0, 0,
                0, 0, scale, 0,
                0, 0, 0, 1
            ]);
        }
        
        // Multiply projection * modelView
        return this.multiplyMatrices(projection, modelView);
    }
    
    multiplyMatrices(a, b) {
        const result = new Float32Array(16);
        for (let i = 0; i < 4; i++) {
            for (let j = 0; j < 4; j++) {
                result[i * 4 + j] = 
                    a[i * 4 + 0] * b[0 * 4 + j] +
                    a[i * 4 + 1] * b[1 * 4 + j] +
                    a[i * 4 + 2] * b[2 * 4 + j] +
                    a[i * 4 + 3] * b[3 * 4 + j];
            }
        }
        return result;
    }
    
    render(pointsData, currentView) {
        const gl = this.gl;
        
        // Update camera for 3D view
        if (currentView === '3d' && this.camera.autoRotate) {
            this.camera.rotationY += this.camera.rotationSpeed;
        }
        
        gl.viewport(0, 0, this.canvas.width, this.canvas.height);
        gl.clearColor(0.0, 0.0, 0.0, 1.0);
        gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);
        
        if (!pointsData || pointsData.length === 0) return;
        
        const numPoints = pointsData.length / 3;
        if (numPoints < 2) return;
        
        // Prepare position and age data
        const positions = new Float32Array(pointsData);
        const ages = new Float32Array(numPoints);
        
        for (let i = 0; i < numPoints; i++) {
            ages[i] = i / (numPoints - 1);
        }
        
        // Create projected positions based on view
        const projectedPositions = new Float32Array(numPoints * 3);
        for (let i = 0; i < numPoints; i++) {
            const x = positions[i * 3];
            const y = positions[i * 3 + 1];
            const z = positions[i * 3 + 2];
            
            switch (currentView) {
                case 'xy':
                    projectedPositions[i * 3] = x;
                    projectedPositions[i * 3 + 1] = y;
                    projectedPositions[i * 3 + 2] = 0;
                    break;
                case 'xz':
                    projectedPositions[i * 3] = x;
                    projectedPositions[i * 3 + 1] = z * 0.5;
                    projectedPositions[i * 3 + 2] = 0;
                    break;
                case 'yz':
                    projectedPositions[i * 3] = y;
                    projectedPositions[i * 3 + 1] = z * 0.5;
                    projectedPositions[i * 3 + 2] = 0;
                    break;
                case '3d':
                default:
                    projectedPositions[i * 3] = x;
                    projectedPositions[i * 3 + 1] = y;
                    projectedPositions[i * 3 + 2] = z;
                    break;
            }
        }
        
        const colors = this.getColors();
        const mvpMatrix = this.createMVPMatrix(currentView);
        
        // Render trails as line strip
        gl.useProgram(this.lineProgram);
        
        gl.bindBuffer(gl.ARRAY_BUFFER, this.positionBuffer);
        gl.bufferData(gl.ARRAY_BUFFER, projectedPositions, gl.DYNAMIC_DRAW);
        gl.enableVertexAttribArray(this.lineAttributes.position);
        gl.vertexAttribPointer(this.lineAttributes.position, 3, gl.FLOAT, false, 0, 0);
        
        gl.bindBuffer(gl.ARRAY_BUFFER, this.ageBuffer);
        gl.bufferData(gl.ARRAY_BUFFER, ages, gl.DYNAMIC_DRAW);
        gl.enableVertexAttribArray(this.lineAttributes.age);
        gl.vertexAttribPointer(this.lineAttributes.age, 1, gl.FLOAT, false, 0, 0);
        
        gl.uniformMatrix4fv(this.lineUniforms.mvpMatrix, false, mvpMatrix);
        gl.uniform3fv(this.lineUniforms.startColor, colors.start);
        gl.uniform3fv(this.lineUniforms.endColor, colors.end);
        
        gl.drawArrays(gl.LINE_STRIP, 0, numPoints);
        
        // Render current point as a larger point
        gl.useProgram(this.pointProgram);
        
        const currentPoint = projectedPositions.slice(-3);
        const currentAge = new Float32Array([1.0]);
        
        gl.bindBuffer(gl.ARRAY_BUFFER, this.positionBuffer);
        gl.bufferData(gl.ARRAY_BUFFER, currentPoint, gl.DYNAMIC_DRAW);
        gl.enableVertexAttribArray(this.pointAttributes.position);
        gl.vertexAttribPointer(this.pointAttributes.position, 3, gl.FLOAT, false, 0, 0);
        
        gl.bindBuffer(gl.ARRAY_BUFFER, this.ageBuffer);
        gl.bufferData(gl.ARRAY_BUFFER, currentAge, gl.DYNAMIC_DRAW);
        gl.enableVertexAttribArray(this.pointAttributes.age);
        gl.vertexAttribPointer(this.pointAttributes.age, 1, gl.FLOAT, false, 0, 0);
        
        gl.uniformMatrix4fv(this.pointUniforms.mvpMatrix, false, mvpMatrix);
        gl.uniform1f(this.pointUniforms.pointSize, currentView === '3d' ? 8.0 : 6.0);
        gl.uniform3fv(this.pointUniforms.startColor, colors.end);
        gl.uniform3fv(this.pointUniforms.endColor, colors.end);
        
        gl.drawArrays(gl.POINTS, 0, 1);
    }
    
    resize(width, height) {
        this.canvas.width = width;
        this.canvas.height = height;
        this.gl.viewport(0, 0, width, height);
    }
    
    setCameraDistance(distance) {
        this.camera.distance = distance;
    }
    
    setCameraRotation(x, y) {
        this.camera.rotationX = x;
        this.camera.rotationY = y;
    }
    
    setAutoRotate(enabled) {
        this.camera.autoRotate = enabled;
    }
    
    setRotationSpeed(speed) {
        this.camera.rotationSpeed = speed;
    }
}