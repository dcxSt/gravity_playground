// let state = []; // the phase space state x1,y1,vx1,vy1,x2,y2,vy2,etc.
// let masses = []; // the masses, must have len = len state / 4
let h = 0.001; // the time step 
let G = 10000;
let canvx = 400;
let canvy = 400;

function setup() {
  frameRate(60);
  createCanvas(canvx, canvy);
  time = 0;
  state = [];
  state.push(new Mass(0,100,900,0,10));
  state.push(new Mass(0,-100,-900,0,12));
  state.push(new Mass(0,0,0,0,10000));
}

function draw() {
  time += h; 
  background(0);
  
  // update
  // update_explicit_euler(state);
  update_rk4(state);
  
  // draw
  state.forEach((mass) => {
    mass.display();
  })
  
}

// takes state: array of Mass objects and updates them by timestep
function update_explicit_euler(state) {
  // compute the derivative
  f_of_y = []; // will be same length as state
  state.forEach((mass, idx) => {
    f_of_y.push(mass.get_gradient(state,idx));
  });
  // update the state
  for (let idx = 0; idx < f_of_y.length; idx ++) {
    let foy = f_of_y[idx];
    let mass = state[idx];
    mass.add_vec(foy); // update the state
  }
}

function divby2(num) {return num/2}
function divby3(num) {return num/3}
// function multby2(num) {return num*2}
function divby6(num) {return num/6}


// takes state: array of Mass objects and updates them by a timestep
function update_rk4(state) {
  // dy_n+1 = y_n + 1/6 h(k1 + k2 + k3 + k4)
  // k1 = f(y_n) // not time-dependent
  // k2 = f(y_n + h*k1 / 2)
  // k3 = f(y_n + h*k2 / 2)
  // k4 = f(y_n + h*k3) 
  
  let k1 = []; // same length as state
  let k2 = []; // same length as state
  let k3 = []; // same length as state
  let k4 = []; // same length as state
  
  // Calculate k1
  state.forEach((mass, idx) => {
    k1.push(mass.get_gradient(state,idx));
  });
  
  // Calculate k2
  let state_copy = deep_copy_state(state);
  // update the state copy
  for (let idx = 0; idx < state.length; idx ++) {
    let mass = state_copy[idx]; // this is a Mass object
    mass.add_vec(k1[idx].map(divby2)); // update the state
  }
  state_copy.forEach((mass, idx) => {
    k2.push(mass.get_gradient(state_copy, idx));
  });
  
  // Calculate k3
  state_copy = deep_copy_state(state);
  // update the state copy
  for (let idx = 0; idx < state.length; idx ++) {
    let mass = state_copy[idx]; // this is a Mass object
    mass.add_vec(k2[idx].map(divby2)); // update the state
  }
  state_copy.forEach((mass, idx) => {
    k3.push(mass.get_gradient(state_copy, idx));
  });
  
  // Calculate k4
  state_copy = deep_copy_state(state);
  // update the state copy
  for (let idx = 0; idx < state.length; idx ++) {
    let mass = state_copy[idx]; // this is a Mass object
    mass.add_vec(k3[idx]); // update the state
  }
  state_copy.forEach((mass, idx) => {
    k4.push(mass.get_gradient(state_copy, idx));
  });
  
  // update the actual state 
  // dy_n+1 = y_n + 1/6 h(k1 + k2 + k3 + k4) 
  for (let idx = 0; idx < state.length; idx ++) {
    let mass = state[idx]; // this is a Mass object
    mass.add_vec(k1[idx].map(divby6));
    mass.add_vec(k2[idx].map(divby3));
    mass.add_vec(k3[idx].map(divby3));
    mass.add_vec(k4[idx].map(divby6))
  }
}

// returns a deep copy of an array of Mass objects
function deep_copy_state() {
  return state.map(deep_copy_mass);
}

// returns deep copy of instance of Mass passed as arg
function deep_copy_mass(mass_obj) {
  return mass_obj.deep_copy();
}

// This computes the acceleration due to potential Hamiltonian 
function get_accel_from(x1,y1,m1,x2,y2,m2) {
  r21 = [x1-x2,y1-y2]; // r pointing from 2 to 1
  rabs = Math.sqrt(r21[0]*r21[0] + r21[1]*r21[1]);
  ax = -G*m2/(rabs*rabs*rabs) * r21[0];
  ay = -G*m2/(rabs*rabs*rabs) * r21[1];
  return [ax,ay];
}

class Mass {
  constructor (x,y,vx,vy,m) {
    this.x = x;
    this.y = y;
    this.vx = vx;
    this.vy = vy;
    this.m = m;
  }
  
  // Explicit Euler; doesn't modify state, just returns value
  get_gradient (state,idx) {
    let ax = 0; // acceleration at this time
    let ay = 0;
    // if it's not the same mass
    state.forEach((item, index) => {
      if (index != idx) {
        let acc = get_accel_from(this.x,this.y,this.m,item.x,item.y,item.m);
        ax += acc[0];
        ay += acc[1];
      }
    });
    // correct formalism dy/dt = f(y)
    let this_f_of_y = [this.vx, this.vy, ax, ay]; // unused, just for reader
    return this_f_of_y; // [dx,dy,dvx,dvy]
    // [this.x,this.y] = [this.x + h*this.vx,this.y + h*this.vy];
    // [this.vx,this.vy] = [this.vx + h*ax,this.vy + h*ay];
  }
  
  add_vec (foy) {
    let [xprime,yprime,vxprime,vyprime] = foy;
    this.x += h * xprime;
    this.y += h * yprime;
    this.vx += h * vxprime;
    this.vy += h * vyprime;
  }
  
  deep_copy () {
    return new Mass(this.x,this.y,this.vx,this.vy,this.m);
  }
  
  display (idx) {
    fill(255 - 40*idx);
    circle(this.x + canvx/2,this.y + canvy/2,min(2*Math.cbrt(this.m),10));
  }
}
