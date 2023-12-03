Node[] nodes;
ArrayList<Edge> edges = new ArrayList<Edge>();
PGraphics g;
String[] names;
float zum = 4.0;

boolean do_text = true;
boolean do_arrows = true;
float node_rad = 14.0;
JSONObject data;
PFont mf;

// "HL", "A"
String col_mode = "HL";

boolean motion = false;
boolean annealing = false;

Gradient cmap;

void setup() {
  size(1024, 1024, P2D);
  noSmooth();
  String[] fl = loadStrings("M.csv");
  String[] coords = loadStrings("coords.csv");
  String[] hl = loadStrings("half_lives.csv");
  cmap = new Gradient("viridis");
  cmap.loadGradient("grad1.png");
  names = loadStrings("names.txt");
  data = loadJSONObject("NUDAT2020.json");
  nodes = new Node[fl.length];
  for(int i = 0; i < fl.length; i++) {
    String[] spl = coords[i].split(" ");
    nodes[i] = new Node();
    nodes[i].A = int(names[i].replaceAll("[^\\d.]", ""));
    if(i < fl.length - 2) {
      if(hl[i] != "inf") {
        nodes[i].hl = float(hl[i]);
      } else {nodes[i].stable = true; nodes[i].hl = 10000000000000000.0;}
    }
    if(true) {
      nodes[i].x = float(spl[0]);
      nodes[i].y = float(spl[1]);
    }
    nodes[i].id = i;
    String[] parts = fl[i].split(",");
    int w = 0;
    for(int j = 0; j < parts.length; j++) {
      if(float(parts[j]) != 0.0 && j != i) {
        w++;
        edges.add(new Edge(i, j));
      }
    }
    if(w > 100) println(w);
  }
  mf = createFont("BarlowCondensed-Medium.ttf", 12);
  
  g = createGraphics(int(8192), int(8192), P2D);
  g.smooth(16);
  
  
  //g = createGraphics(1024, 1024, P2D);
  //g.noSmooth();
}
void draw() {
  g.beginDraw();
  g.background(0);
  g.stroke(255);
  //g.blendMode(ADD);
  g.translate(2048*zum/4.0, 2048*zum/4.0);
  g.fill(0);
  g.textAlign(CENTER, CENTER);
  g.textSize(12);
  g.textFont(mf);
  g.strokeWeight(node_rad*2);
  for(Node n : nodes) {
    n.disp();
    if(motion) n.move();
  }
  g.fill(0);
  g.strokeWeight(1);
  for(Edge e : edges) {
    if(e.i < nodes.length-2 && e.j < nodes.length-2) {
      e.disp();
      if(motion) e.force();
    }
  }
  if(motion)
  if(annealing)
  for(int i = 0; i < 100000; i++) {
    repulse((int)random(nodes.length), (int)random(nodes.length));
  }
  else
  for(int i = 0; i < nodes.length; i++) {
    for(int j = 0; j < i; j++) {
      repulse(i, j);
    }
  }
  if(motion)
  for(int i = 0; i < nodes.length; i++) {
    mousePulse(i);
    inPulse(i);
  }
  g.blendMode(BLEND);
  g.endDraw();
  image(g, 0, 0, height, height);
  //if(mousePressed) strokeWeight(2); else strokeWeight(1);
}

void keyPressed() {
  if(key == 'p') {
    g.save("portrait.png");
    print("saved");
  }
  if(key == ' ') {
    String[] coords = new String[nodes.length];
    for(int i = 0; i < nodes.length; i++) {
      coords[i] = nodes[i].x + " " + nodes[i].y;
    }
    saveStrings("coords.csv", coords);
  }
}
void mousePulse(int i) {
  if(mousePressed) {
    float dx = nodes[i].x - mouseX;
    float dy = nodes[i].y - mouseY;
    float r = sqrt(dx*dx + dy*dy);
    float f = -10.0/(pow(r, 1.7)+10.0);
    float fx = dx * f;
    float fy = dy * f;
    nodes[i].xv -= fx;
    nodes[i].yv -= fy;
  }
}

void inPulse(int i) {
  float dx = nodes[i].x - 512;
  float dy = nodes[i].y - 512;
  float r = sqrt(dx*dx + dy*dy);
  float f = 0.005/(pow(r, 0.5)+10.0);
  float fx = dx * f;
  float fy = dy * f;
  nodes[i].xv -= fx;
  nodes[i].yv -= fy;
}

void repulse(int i, int j) {
  float dx = nodes[i].x - nodes[j].x;
    float dy = nodes[i].y - nodes[j].y;
    float r = sqrt(dx*dx + dy*dy);
    float f = -(annealing?40.8:0.8)/(pow(r, 2.5)+4.0);
    if(r < 2.5*node_rad) r *= 6.0*abs(r - 2.5*node_rad)/node_rad + 1.0;
    float fx = dx * f;
    float fy = dy * f;
    nodes[i].xv -= fx;
    nodes[i].yv -= fy;
    nodes[j].xv += fx;
    nodes[j].yv += fy;
}

class Edge {
  int i;
  int j;
  color c;
  Edge(int i, int j) {
    this.i = j;
    this.j = i;
    c = color(100, 30, 255);
  }
  void disp() {
    if(nodes[i].A == nodes[j].A) c = color(255);
    else if(nodes[i].A == nodes[j].A - 4) c = color(40, 255, 180);
    
    g.fill(c);
    g.stroke(c);
    if(do_arrows)
    drawArrow(nodes[j].x*zum, nodes[j].y*zum, nodes[i].x*zum, nodes[i].y*zum, 10, 5, node_rad);
    else 
    g.line(nodes[j].x*zum, nodes[j].y*zum, nodes[i].x*zum, nodes[i].y*zum);
  }
  void force() {
    float dx = nodes[i].x - nodes[j].x;
    float dy = nodes[i].y - nodes[j].y;
    float r = sqrt(dx*dx + dy*dy);
    float f = 0.001*(r - 10.0);
    float fx = dx * f;
    float fy = dy * f;
    nodes[i].xv -= fx;
    nodes[i].yv -= fy;
    nodes[j].xv += fx;
    nodes[j].yv += fy;
  }
}

color HL_to_color(float f) {
  float skinch = log(f)/log(10);
  return cmap.getGrad(skinch, 10.0, -8.0);
}

float luma(color c) {
  return 0.21*red(c) + 0.71*green(c) + 0.07*blue(c);
}

class Node {
  float x = 0;
  float y = 0;
  float xv = 0;
  float yv = 0;
  int id = 0;
  float hl = 0.0;
  int A = 0;
  boolean stable = false;
  Node() {
    x = width/2 + random(-10, 10);
    y = height/2 + random(-10, 10);
  }
  void disp() {
    g.fill(0);
    if(hl == 0.0) hl = 10000000000000000000.0;
    color c = color(255);
    switch(col_mode) {
      case "HL":
        if(!stable) c = HL_to_color(hl);
        break;
      case "A":
        c = cmap.getGrad((float)A, 0.0, 294.0);
        break;
    }
    g.stroke(c);
    g.point(x*zum, y*zum);
    if(do_text) {
     if(luma(c) > 255/2.f) g.fill(0); else g.fill(255);
      g.text(names[id], x*zum, y*zum - 2);
    }
  }
  void move() {
    x += xv*0.4;
    y += yv*0.4;
    xv *= 0.97;
    yv *= 0.97;
  }
}

void drawArrow(float x00, float y00, float x11, float y11, float headLength, float headWidth, float offrad) {
  PVector d = new PVector(x11 - x00, y11 - y00);
  d.normalize();
  d.mult(offrad);
  float x0 = x00 + d.x;
  float y0 = y00 + d.y;
  float x1 = x11 - d.x;
  float y1 = y11 - d.y;
  g.line(x0, y0, x1, y1);
  PVector v = new PVector(x1 - x0, y1 - y0);
  v.normalize();
  PVector ortho = new PVector(-v.y, v.x);
  v.mult(headLength);
  ortho.mult(headWidth);
  g.noStroke();
  g.triangle(x1, y1, x1 - v.x + ortho.x, y1 - v.y + ortho.y, x1 - v.x - ortho.x, y1 - v.y - ortho.y);
  g.stroke(255);
}
