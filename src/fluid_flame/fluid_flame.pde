float videoScale;

// Number of columns and rows in our system
int cols, rows;

int N = 64;

int SIZE = (N+2)*(N+2);

//v is velocities
float[] u = new float[SIZE];
float[] v = new float[SIZE];
float[] u_prev = new float[SIZE];
float[] v_prev = new float[SIZE];

float[] dens = new float[SIZE];
float[] dens_prev = new float[SIZE];

//Entry is true if the grid contains an object
boolean[] object_space_grid = new boolean[SIZE];

float force = 5.0f;
float source = 100.0f;
float dt = 0.4f;
float visc = 0.0f;
float diff = 0.0f;

//For ui
boolean mousePressed = false;
boolean mLeftPressed = false;
boolean mRightPressed = true;
int mx, my; //Current mouse position
int omx, omy; //Old mouse position
boolean drawVel = false;
boolean advectNormalStep = true;
boolean object_space = false;
boolean fancyColours = true;

//Helper function to get an element from a 1D array as if it were a 2D array
int IX(int i, int j)
{
    return (i + ((N+2) * j));
}

void clearData()
{
    int i;
    int sz = (N+2)*(N+2);

    for(i=0;i<sz; i++)
    {
        u[i] = 0.0f;
        v[i] = 0.0f;
        u_prev[i] = 0.0f;
        v_prev[i] = 0.0f;
        dens[i] = 0.0f;
        dens_prev[i] = 0.0f;
        object_space_grid[i] = false;
    }

    //Set initial object spaces
/*
  object_space_grid[IX(10,10)] = true;
  object_space_grid[IX(11,10)] = true;
  object_space_grid[IX(12,10)] = true;
  object_space_grid[IX(13,10)] = true;

  object_space_grid[IX(10,11)] = true;
  object_space_grid[IX(11,11)] = true;
  object_space_grid[IX(12,11)] = true;
  object_space_grid[IX(13,11)] = true;

  object_space_grid[IX(10,12)] = true;
  object_space_grid[IX(11,12)] = true;
  object_space_grid[IX(12,12)] = true;
  object_space_grid[IX(13,12)] = true;

  object_space_grid[IX(10,13)] = true;
  object_space_grid[IX(11,13)] = true;
  object_space_grid[IX(12,13)] = true;
  object_space_grid[IX(13,13)] = true;
*/

    int j;
    for(i=N/4;i<=(3*N/4);i++)
    {
        for(j=N/4;j<=(3*N/4);j++)
        {
            object_space_grid[IX(i,j)] = true;
        }
    }

}

void setup() {
    size(512,512, P3D);
    background(0);
//    frameRate(30);
//    lights();
    //   ortho(0.0f, 1.0f, 1.0f, 0.0f, -10.0f, 10.0f);
//    ortho(0,1, 0,1, -10,20);

    // Initialize columns and rows
//  cols = width/videoScale;
//  rows = height/videoScale;
    videoScale = height / (N+2);
    rows = N+2;
    cols = N+2;
    clearData();
}

void drawRectangles()
{
    // Begin loop for columns
    for (int i = 0; i < cols; i++) {
        // Begin loop for rows
        for (int j = 0; j < rows; j++) {

            // Scaling up to draw a rectangle at (x,y)
            float x = i*videoScale;
            float y = j*videoScale;
            fill(255);
            //fill(243,32,32);
            stroke(0);
            // For every column and row, a rectangle is drawn at an (x,y) location scaled and sized by videoScale.
            rect(x,y,videoScale,videoScale);
        }
    }
}

void drawSquare(float x, float y, float side)
{
    rect(x*width,y*height,side*width,side*height);
}

void drawObjects()
{
    stroke(0);
    fill(123);

    float x, y, h;
    h = 1.0f / (float)N;
    for(int i = 0; i < N+2 ; i++)
    {
        x = (i-0.5f)*h;
        for(int j = 0 ; j < N+2 ;j++)
        {
            y = (j-0.5f)*h;
            if(object_space_grid[IX(i,j)])
            {
                drawSquare(x,y,h);
            }
        }
    }
}

/**
 * Sets boundary for diffusion. It is bound vertically and horizontally in a box.
 **/
void setBnd(int n, int b, float[] x)
{
    int i,j;

    if(object_space)
    {
        //Bound along internal boundaries
        for(i=1;i<=n;i++)
        {
            for(j=1;j<=n;j++)
            {
                if(object_space_grid[IX(i,j)])
                {
                    //x[IX(i,j)] = (((x[IX(i - 1, j)] + x[IX(i + 1,j)] + x[IX(i, j - 1)] + x[IX(i,j + 1)])) / 4);

                    x[IX(i,j)] = x[IX(i-1,j)];
                    x[IX(i,j)] -= x[IX(i,j-1)];
                    x[IX(i,j)] += x[IX(i+1,j)];
                    x[IX(i,j)] -= x[IX(i+1,j+1)];

                    x[IX(i,j)] *= 0.25;
                }
            }
        }
    }


    //Bound along edges, will wrap around level if the flag is set
    for(i=0 ; i<=n ; i++)
    {
        x[IX(0,i)] = (b==1 ? -x[IX(1,i)] : x[IX(1,i)]);
        x[IX(n+1,i)] = (b==1 ? -x[IX(n,i)] : x[IX(n,i)]);
        x[IX(i,0)] = (b==2 ? -x[IX(i,1)] : x[IX(i,1)]);
        x[IX(i,n+1)] = (b==2 ? -x[IX(i,n)] : x[IX(i,n)]);

        if(! advectNormalStep)
        {
            //Reflect in on itself
            float temp = x[IX(0,i)];
            x[IX(0,i)] = x[IX(n+1,i)];
            x[IX(n+1,i)] = temp;

            temp = x[IX(i,0)];
            x[IX(i,0)] = x[IX(i,n+1)];
            x[IX(i,n+1)] = temp;
        }
/*
  x[IX(0,i)] = (b==1 ? -x[IX(1,i)] : x[IX(1,i)]);
  x[IX(n+1,i)] = (b==1 ? -x[IX(n,i)] : x[IX(n,i)]);
  x[IX(i,0)] = (b==2 ? -x[IX(i,1)] : x[IX(i,1)]);
  x[IX(i,n+1)] = (b==2 ? -x[IX(i,n)] : x[IX(i,n)]);
*/
    }
    x[IX(0,0)] = 0.5f * (x[IX(1,0)] + x[IX(0,1)]);
    x[IX(0,n+1)] = 0.5f * (x[IX(1,n+1)] + x[IX(0,n)]);
    x[IX(n+1,0)] = 0.5f * (x[IX(n,0)] + x[IX(n+1,1)]);
    x[IX(n+1,n+1)] = 0.5f * (x[IX(n,n+1)] + x[IX(n+1 ,n)]);
}
/*
  void setBnd(int n, int b, float[] x)
  {
  int i;

  for(i=0 ; i<=n ; i++)
  {
  x[IX(0,i)] = (b==1 ? -x[IX(1,i)] : x[IX(1,i)]);
  x[IX(n+1,i)] = (b==1 ? -x[IX(n,i)] : x[IX(n,i)]);
  x[IX(i,0)] = (b==2 ? -x[IX(i,1)] : x[IX(i,1)]);
  x[IX(i,n+1)] = (b==2 ? -x[IX(i,n)] : x[IX(i,n)]);
  }
  x[IX(0,0)] = 0.5f * (x[IX(1,0)] + x[IX(0,1)]);
  x[IX(0,n+1)] = 0.5f * (x[IX(1,n+1)] + x[IX(0,n)]);
  x[IX(n+1,0)] = 0.5f * (x[IX(n,0)] + x[IX(n+1,1)]);
  x[IX(n+1,n+1)] = 0.5f * (x[IX(n,n+1)] + x[IX(n+1 ,n)]);
  }
*/
/**
 * 1st step:
 * Add a source held in s to the density held in x.
 * dt is the timestep and n+2 is the row and col size.
 */
void addSource(int n, float[] x, float[] s, float dt)
{
    int i;
    int sz = (n+2)*(n+2);

    for(i=0;i<sz;i++)
    {
        x[i] += dt*s[i];
    }
}

/**
 * 2nd Step:
 * Diffusion at rate diff. Each cell will exchange density with direct neighbours.
 * Uses Gauss-Seidel relaxation.
 */
void diffuse(int n, int b,  float[] x, float[] x0, float diff, float dt)
{
    int i, j, k;
    float a = dt * diff*n*n;

    if(advectNormalStep)
    {
        for(k=0 ; k<20; k++)
        {
            for(i=1;i<=n;i++)
            {
                for(j=1;j<=n;j++)
                {


                    x[IX(i,j)] = (x0[IX(i,j)] +
                                  a * (x[IX(i - 1, j)] + x[IX(i + 1,j)] + x[IX(i, j - 1)] + x[IX(i,j + 1)])) / (1 + (4 * a));

                }
            }
            setBnd(n,b,x);
        }
    }
    else
    {
        for(k=0 ; k<20; k++)
        {
            for(i=0;i<=n+1;i++)
            {
                for(j=0;j<=n+1;j++)
                {
                    //Need to check for boundaries
                    if(i==0)
                    {
                        if(j==0)
                        {
                            x[IX(i,j)] = (x0[IX(i,j)] +
                                          a * (x[IX(n+1, j)] + x[IX(i+1,j)] + x[IX(i, n+1)] + x[IX(i,j + 1)])) / (1 + (4 * a));
                        }
                        else if(j==n+1)
                        {
                            x[IX(i,j)] = (x0[IX(i,j)] +
                                          a * (x[IX(n+1, j)] + x[IX(i + 1,j)] + x[IX(i, j - 1)] + x[IX(i,0)])) / (1 + (4 * a));
                        }
                        else
                        {
                            x[IX(i,j)] = (x0[IX(i,j)] +
                                          a * (x[IX(n+1, j)] + x[IX(i + 1,j)] + x[IX(i, j - 1)] + x[IX(i,j + 1)])) / (1 + (4 * a));
                        }
                    }
                    else if(i==n+1)
                    {
                        if(j==0)
                        {
                            x[IX(i,j)] = (x0[IX(i,j)] +
                                          a * (x[IX(i-1, j)] + x[IX(0,j)] + x[IX(i, n+1)] + x[IX(i,j + 1)])) / (1 + (4 * a));
                        }
                        else if(j==n+1)
                        {
                            x[IX(i,j)] = (x0[IX(i,j)] +
                                          a * (x[IX(i-1, j)] + x[IX(0,j)] + x[IX(i, j-1)] + x[IX(i,0)])) / (1 + (4 * a));
                        }
                        else
                        {
                            x[IX(i,j)] = (x0[IX(i,j)] +
                                          a * (x[IX(i-1, j)] + x[IX(0,j)] + x[IX(i, j - 1)] + x[IX(i,j + 1)])) / (1 + (4 * a));
                        }
                    }
                    else
                    {
                        if(j==0)
                        {
                            x[IX(i,j)] = (x0[IX(i,j)] +
                                          a * (x[IX(i - 1, j)] + x[IX(i + 1,j)] + x[IX(i, n+1)] + x[IX(i,j + 1)])) / (1 + (4 * a));
                        }
                        else if(j==n+1)
                        {
                            x[IX(i,j)] = (x0[IX(i,j)] +
                                          a * (x[IX(i - 1, j)] + x[IX(i + 1,j)] + x[IX(i, j - 1)] + x[IX(i,0)])) / (1 + (4 * a));
                        }
                        else
                        {
                            x[IX(i,j)] = (x0[IX(i,j)] +
                                          a * (x[IX(i - 1, j)] + x[IX(i + 1,j)] + x[IX(i, j - 1)] + x[IX(i,j + 1)])) / (1 + (4 * a));
                        }
                    }
                }
            }
            setBnd(n,b,x);
        }
    }
}

void advectBoundaryStep(int n,
                        float[] d,
                        float[] d0,
                        float[] u,
                        float[] v,
                        float dt,
                        float dt0,
                        int i,
                        int j)
{
    int i0, i1, j0, j1;
    float x, y, s0, t0, s1, t1;

    x = i - dt0 * u[IX(i,j)];
    y = j - dt0 * v[IX(i,j)];

    if(x < 0.5f)
    {
        x = n + 0.5f;
    }
    if(x > (n+0.5f))
    {
        x = 0.5f;
    }
    i0 = (int)x;
    i1 = i0 + 1;

    if(y < 0.5f)
    {
        y = n + 0.5f;
    }
    if(y > (n+0.5f))
    {
        y = 0.5f;
    }
    j0 = (int)y;
    j1 = j0 + 1;

    s1 = x - i0;
    s0 = 1 - s1;

    t1 = y - j0;
    t0 = 1 - t1;

    d[IX(i,j)] = s0 * (t0 * d0[IX(i0,j0)] + t1 * d0[IX(i0,j1)]) +
        s1 * (t0 * d0[IX(i1,j0)] + t1 * d0[IX(i1,j1)]);
}

void advectNormalStep(int n,
                      float[] d,
                      float[] d0,
                      float[] u,
                      float[] v,
                      float dt,
                      float dt0,
                      int i,
                      int j)
{
    int i0, i1, j0, j1;
    float x, y, s0, t0, s1, t1;

    x = i - dt0 * u[IX(i,j)];
    y = j - dt0 * v[IX(i,j)];

    if(x < 0.5f)
    {
        x = 0.5f;
    }
    if(x > (n+0.5f))
    {
        x = n + 0.5f;
    }
    i0 = (int)x;
    i1 = i0 + 1;

    if(y < 0.5f)
    {
        y = 0.5f;
    }
    if(y > (n+0.5f))
    {
        y = n + 0.5f;
    }
    j0 = (int)y;
    j1 = j0 + 1;

    s1 = x - i0;
    s0 = 1 - s1;

    t1 = y - j0;
    t0 = 1 - t1;

    d[IX(i,j)] = s0 * (t0 * d0[IX(i0,j0)] + t1 * d0[IX(i0,j1)]) +
        s1 * (t0 * d0[IX(i1,j0)] + t1 * d0[IX(i1,j1)]);
}

/**
 * 3rd Step :
 * Calculating advections. This ensures that the density follows a given velocity field.
 **/
void advect(int n,
            int b,
            float[] d,
            float[] d0,
            float[] u,
            float[] v,
            float dt)
{
    int i, j, i0, j0, i1, j1;
    float x, y, s0, t0, s1, t1, dt0;

    dt0 = dt*n;
    if(advectNormalStep)
    {
        for(i=1 ; i<=n ; i++)
        {
            for(j=1 ; j<=n ; j++)
            {
                advectNormalStep(n, d, d0, u, v, dt, dt0, i, j);
            }
        }
    }
    else
    {
        for(i=1 ; i<=n ; i++)
        {
            for(j=1 ; j<=n ; j++)
            {
                advectBoundaryStep(n, d, d0, u, v, dt, dt0, i, j);
            }
        }
    }
    setBnd(n, b, d);
}
/*
  void swap(float[] x, float[] y)
  {
  float[] temp = x;
  x = y;
  y = temp;
  }
*/

/**
 * 1 step of the density solver.
 */
void densStep(int n, float[] x, float[] x0, float[] u, float[] v, float diff, float dt)
{
    addSource(n, x, x0, dt);
    float[] temp = x0;
    x0 = x;
    x = temp;
//   swap(x0, x);
    diffuse(n, 0, x, x0, diff, dt);
//   swap(x0, x);
    temp = x0;
    x0 = x;
    x = temp;
    advect(n, 0, x, x0, u, v, dt);
}

/**
 * 1 step of the velocity solver.
 */
void velStep( int n,
              float[] u,
              float[] v,
              float[] u0,
              float[] v0,
              float visc,
              float dt)
{
    addSource(n, u, u0, dt);
    addSource(n, v, v0, dt);

//    swap(u0, u);
    float[] temp = u0;
    u0 = u;
    u = temp;

    diffuse(n, 1, u, u0, visc, dt);

//    swap(v0, v);
    temp = v0;
    v0 = v;
    v = temp;

    diffuse(n, 2, v, v0, visc, dt);
    project(n, u, v, u0, v0);

//    swap(u0, u);
    temp = u0;
    u0 = u;
    u = temp;

//    swap(v0, v);
    temp = v0;
    v0 = v;
    v = temp;

    advect(n, 1, u, u0, u0, v0, dt);
    advect(n, 2, v, v0, u0, v0, dt);
    project(n, u, v, u0, v0);
}

void project(int n,
             float[] u,
             float[] v,
             float[] p,
             float[] div)
{
    int i,j,k;
    float h = 1.0/n;
    for(i=1; i<=n ;i++)
    {
        for(j=1;j<=n;j++)
        {
            div[IX(i,j)] = -0.5f*h*(u[IX(i+1,j)] - u[IX(i-1,j)] +
                                    v[IX(i,j+1)] - v[IX(i,j-1)]);
            p[IX(i,j)] = 0;
        }
    }
    setBnd(n, 0, div);
    setBnd(n, 0, p);

    for(k=0;k<20;k++)
    {
        for(i=1;i<=n;i++)
        {
            for(j=1;j<=n;j++)
            {
                p[IX(i,j)] = (div[IX(i,j)]+p[IX(i-1,j)]+p[IX(i+1,j)]+
                              p[IX(i,j-1)]+p[IX(i,j+1)]) / 4;
            }
        }
        setBnd(n, 0, p);
    }

    for(i=1;i<=n;i++)
    {
        for(j=1;j<=n;j++)
        {
            u[IX(i,j)] -= 0.5f*(p[IX(i+1,j)]-p[IX(i-1,j)])/h;
            v[IX(i,j)] -= 0.5f*(p[IX(i,j+1)]-p[IX(i,j-1)])/h;
        }
    }
    setBnd(n, 1, u);
    setBnd(n, 2, v);
}

void drawPointOnScreen(float x, float y)
{
    vertex(x * width, y * height);
}

void colourFancy(float colour, float threshold)
{
    if(colour < threshold)
    {
	stroke(colour, 0, 0);
	fill(colour, 0, 0);
    }
    else
    {
	stroke(colour, colour/2, 0);
	fill(colour, colour/2, 0);
    }
}

void drawDensity()
{
    int i,j;
    float x, y, h, d00, d01, d10, d11;

    h = 1.0f / (float)N;

    for(i=0;i<=N;i++)
    {
        x = (i-0.5f)*h;

        for(j=0;j<=N;j++)
        {
            y = (j-0.5f)*h;

            d00 = dens[IX(i,j)];
            d01 = dens[IX(i,j+1)];
            d10 = dens[IX(i+1,j)];
            d11 = dens[IX(i+1,j+1)];

            d00 *= 255;
            d10 *= 255;
            d11 *= 255;
            d01 *= 255;

            beginShape(QUADS);


            float threshold = 200.0f;
            if(fancyColours)
            {
		colourFancy(d00, threshold);
            }
            else
            {
                stroke(d00, d00, d00);
                fill(d00, d00, d00);
            }
            drawPointOnScreen(x, y);


            if(fancyColours)
            {
		colourFancy(d10, threshold);
            }
            else
            {
                stroke(d10, d10, d10);
                fill(d10, d10, d10);
            }
            drawPointOnScreen(x+h, y);



            if(fancyColours)
            {
		colourFancy(d11, threshold);
            }
            else
            {
                stroke(d11,d11,d11);
                fill(d11, d11, d11);
            }
            drawPointOnScreen(x+h, y+h);

            if(fancyColours)
            {
		colourFancy(d01, threshold);
            }
            else
            {
                stroke(d01,d01,d01);
                fill(d01, d01, d01);
            }
            drawPointOnScreen(x, y+h);

            endShape();
        }
    }
}

void drawLineOnScreen(float x1, float y1, float x2, float y2)
{
    line(x1*width, y1*height,
         x2*width, y2*height);
}

void drawVelocity()
{
    int i,j;
    float x, y, h;

    h = 1.0f / (float)N;

    for(i=1;i<=N;i++)
    {
        x = (i-0.5f)*h;

        for(j=1;j<=N;j++)
        {
            y = (j-0.5f)*h;

            stroke(255);
            fill(255);

            drawLineOnScreen(x, y, x + u[IX(i,j)], y + v[IX(i,j)]);
        }
    }
}

void getForcesFromUI(float[] d, float[] u, float[] v)
{
    int i, j;
    int size = (N+2)*(N+2);

    for(i=0;i<size;i++)
    {
        u[i] = 0.0f;
        v[i] = 0.0f;
        d[i] = 0.0f;
    }

    if(! mousePressed)
    {
        return;
    }

    i = (int)((mx / (float)width)*N+1);
    j = (int)((my / (float)height)*N+1);

    if( (i<1) || (i>N)
        || (j<1) || (j>N))
    {
        return;
    }

    if(mLeftPressed)
    {
        u[IX(i,j)] = force * (mx-omx);
        v[IX(i,j)] = -force * (omy-my);
    }

    if(mRightPressed)
    {
        d[IX(i,j)] = source; //Set density to initial value
    }

    omx = mx;
    omy = my;
}

void mousePressed() {
    mousePressed = true;

    omx = mouseX;
    mx = mouseX;
    omy = mouseY;
    my = mouseY;

    if(mouseButton == LEFT)
    {
        mLeftPressed = true;
    }
    else if (mouseButton == RIGHT)
    {
        mRightPressed = true;
    }
}

void mouseReleased()
{
    mousePressed = false;

    if(mouseButton == LEFT)
    {
        mLeftPressed = false;
    }
    else if (mouseButton == RIGHT)
    {
        mRightPressed = false;
    }
}

void keyPressed()
{
    if (key == 'c' || key == 'C')
    {
        clearData();
    }
    if( key == 'v' || key == 'V')
    {
        println("v pressed");
        drawVel = ! drawVel;
    }
    if( key == 'a' || key == 'A')
    {
        advectNormalStep = ! advectNormalStep;
    }
    if( key == 'o' || key == 'O')
    {
        object_space = ! object_space;
    }
    if( key == 'f' || key =='F')
    {
        fancyColours = ! fancyColours;
    }
}

void draw()
{
    mx = mouseX;
    my = mouseY;

/*
  if(mouseButton == LEFT)
  {
  mLeftPressed = true;
  }
  else if (mouseButton == RIGHT)
  {
  mRightPressed = true;
  }
*/

    getForcesFromUI(dens_prev, u_prev, v_prev);
    velStep(N, u, v, u_prev, v_prev, visc, dt);
    densStep(N, dens, dens_prev, u, v, diff, dt);

//    drawRectangles();

    background(0);

    if(drawVel)
    {
        drawVelocity();
    }
    else
    {
        drawDensity();
    }

    if(object_space)
    {
        drawObjects();
    }
}


