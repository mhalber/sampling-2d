// compile: gcc -g3 sampling.c -o sampling -lglfw3 -lopengl32 -lglew32 -lnanovg

#include <stdio.h>
#include <math.h>
#include <stdint.h>
#include <assert.h>
#include <GL/glew.h>
#define GLFW_INCLUDE_GLEXT
#include <GLFW/glfw3.h>
#include "nanovg/nanovg.h"
#define NANOVG_GL3_IMPLEMENTATION
#include "nanovg/nanovg_gl.h"

static const float PADDING = 25.0f;

float randf()
{
  return (float)(rand()) / (float)(RAND_MAX);
}

float randf_range( float min, float max )
{
  return min + (float)(rand()) / (float)(RAND_MAX/(max-min));
}

float
rand_sign()
{
  return ((rand() % 2) * 2 - 1);
}

//=================================================================================================
typedef struct vec2
{
  float x;
  float y;
} vec2_t;

vec2_t vec2_add( const vec2_t p, const vec2_t q )
{
  vec2_t w = (vec2_t){ .x = p.x + q.x, .y = p.y + q.y };
  return w;
}

vec2_t vec2_sub( const vec2_t p, const vec2_t q )
{
  vec2_t w = (vec2_t){ .x = p.x - q.x, .y = p.y - q.y };
  return w;
}

float vec2_norm( const vec2_t p )
{
  return sqrtf(p.x * p.x + p.y * p.y);
}

float vec2_dist( const vec2_t p, const vec2_t q )
{
  vec2_t w = { .x = p.x - q.x, .y = p.y - q.y };
  return sqrtf(w.x * w.x + w.y * w.y);
}

//==================================================================================================
vec2_t* 
random_sampling( size_t n_pts, vec2_t domain_size )
{
  vec2_t* pts = malloc( n_pts * sizeof(vec2_t) );
  for( size_t i = 0; i < n_pts; ++i )
  {
    pts[i] = (vec2_t){ .x = (randf() * domain_size.x),
                       .y = (randf() * domain_size.x) };
  }
  for( size_t i = 0; i < n_pts; ++i )
  {
    pts[i].x += PADDING;
    pts[i].y += PADDING;
  }
  return pts;
}

vec2_t* 
stratified_sampling( vec2_t domain_size, float cell_size, int *n_pts )
{
  int res_x = domain_size.x / cell_size;
  int res_y = domain_size.y / cell_size;
  *n_pts = res_x*res_y;
  vec2_t* pts = malloc((*n_pts)*sizeof(vec2_t));
  
  int idx = 0;
  for( size_t y = 0; y < res_y; y++ )
  {
    for( size_t x = 0; x < res_x; x++ )
    {
      float offset = randf();
      pts[idx++] = (vec2_t){ .x = (x + offset) * cell_size,
                             .y = (y + offset) * cell_size };
    }
  }
  for( size_t i = 0; i < *n_pts; ++i )
  {
    pts[i].x += PADDING;
    pts[i].y += PADDING;
  }
  return pts;
}

vec2_t*
generate_circular_samples_range( const vec2_t c, const float r1, const float r2, const size_t k )
{
  int n = 0;
  vec2_t* pts = malloc( k * sizeof(vec2_t) );
  while( n != k )
  {
    vec2_t offset = { .x = (2 * randf() - 1) * r2,
                      .y = (2 * randf() - 1) * r2 };
    vec2_t candidate = vec2_add( c, offset );
    float dist = vec2_dist( candidate, c );
    if( dist >= r1 && dist < r2 )
    {
      pts[n++] = candidate;
    }
  }
  return pts;
}

int
pt_to_idx( vec2_t pt, int res_y, float cell_size )
{
  int x = pt.x / cell_size;
  int y = pt.y / cell_size;
  return y * res_y + x;
}

vec2_t*
poisson_disk_sampling( vec2_t domain_size, float radius, int k, int* n_pts )
{
  // step 0. prepare domain
  float cell_size = radius / sqrt(2);
  int res_x = domain_size.x / cell_size;
  int res_y = domain_size.x / cell_size;
  int max_n_pts = res_x * res_y;
  vec2_t *pts = malloc( max_n_pts * sizeof(vec2_t) );
  int *indices = malloc( res_x * res_y * sizeof(int) );
  int *active_list = malloc( max_n_pts * sizeof(int) );
  for( size_t i = 0; i < res_x * res_y; ++i )
  {
    indices[i] = -1;
  }

  // step 1. generate random sample
  int pt_idx = 0;
  int al_idx = 0;
  active_list[al_idx++] = pt_idx;
  pts[pt_idx] = (vec2_t){ .x = (randf() * domain_size.x),
                          .y = (randf() * domain_size.y) };
  int idx = pt_to_idx( pts[pt_idx], res_y, cell_size );
  indices[idx] = pt_idx;
  pt_idx++;

  // step 2. to do;
  int counter = 0;
  while( al_idx > 0 )
  {
    int found = 0;
    vec2_t prev_pt = pts[ active_list[al_idx-1] ];
    vec2_t* candidates = generate_circular_samples_range( prev_pt, radius, 2*radius, k );
    for( int i = 0; i < k; ++i )
    {
      // check if candidate is within domain
      if( candidates[i].x < 0 || candidates[i].x >= domain_size.x || 
          candidates[i].y < 0 || candidates[i].y >= domain_size.y ) { continue; }
      
      // check if candidate is using free bin
      int idx_x = candidates[i].x / cell_size;
      int idx_y = candidates[i].y / cell_size;
      if( idx_x >= res_x || idx_y >= res_y ) { continue; }
      int cand_idx = idx_y * res_y + idx_x;
      if( indices[cand_idx] >= 0 ) { continue; }


      // do neighbor search to see if no samples closer that 'radius'
      int invalid = 0;
      for( int oy = -3; oy <= 3; ++oy )
      {
        int iy = idx_y + oy;
        if( iy < 0 || iy >= res_y ) { continue; }
        for( int ox = -3; ox <= 3; ++ox )
        {
          int ix = idx_x + ox;
          if( ix < 0 || ix >= res_x ) { continue; }
          int neigh_idx = iy * res_y + ix;
          assert( neigh_idx >= 0 && neigh_idx < res_x*res_y );
          if( indices[neigh_idx] >= 0 )
          {
            vec2_t neigh = pts[indices[neigh_idx]];
            float dist = vec2_dist( candidates[i], neigh );
            if( dist < radius ) { invalid = 1; break; }
          }
        }
        if( invalid ) { break; }
      }

      if( !invalid )
      {
        pts[pt_idx] = candidates[i];
        active_list[al_idx++] = pt_idx;
        indices[cand_idx] = pt_idx;
        pt_idx++;
        found = 1;
        break;
      }
    }
    if( !found ) { al_idx--; }
    free(candidates);
  }

  // step3. finish up
  free( indices );
  free( active_list );
  *n_pts = pt_idx;
  for( size_t i = 0; i < *n_pts; ++i )
  {
    pts[i].x += PADDING;
    pts[i].y += PADDING;
  }
  return pts;
}

void test_poisson_samples( vec2_t* pts, size_t n_pts, float radius )
{
  float min_dist = 1e9;
  for( size_t i = 0; i < n_pts; ++i )
  {
    vec2_t p = pts[i];
    for( size_t j = 0; j < n_pts; ++j )
    {
      if( i == j ) { continue; }
      vec2_t q = pts[j];
      float dist = vec2_dist(p, q);
      if( dist < min_dist ) { min_dist = dist; }
      if( dist < radius ) 
      { 
        printf("Test Failed!\n"); 
        printf("(%f, %f) (%f, %f) -> %f | %f\n", p.x, p.y, q.x, q.y, dist, radius );
        return; 
      }
    }
  }
  printf("Test passed! Min Dist: %f\n", min_dist );
}

//=================================================================================================

void 
draw_grid( NVGcontext* vg, const vec2_t domain_size, const float cell_size )
{
  int res_x = domain_size.x / cell_size;
  int res_y = domain_size.y / cell_size;
  nvgBeginPath( vg );
  for( size_t i = 0; i < res_x+1; i++ )
  {
    float ox = i * cell_size;
    vec2_t p1 = (vec2_t){ .x = ox + PADDING, .y = PADDING };
    vec2_t p2 = (vec2_t){ .x = ox + PADDING, .y = domain_size.y + PADDING };
    nvgMoveTo( vg, p1.x, p1.y );
    nvgLineTo( vg, p2.x, p2.y );
  }
  for( size_t i = 0; i < res_y+1; ++i )
  {
    float oy = i * cell_size;
    vec2_t p1 = (vec2_t){ .x = PADDING, .y = oy + PADDING };
    vec2_t p2 = (vec2_t){ .x = domain_size.x + PADDING, .y = oy + PADDING };
    nvgMoveTo( vg, p1.x, p1.y );
    nvgLineTo( vg, p2.x, p2.y );
  }

  nvgStrokeColor( vg, nvgRGBA( 128, 128, 128, 255 ) );
  nvgStrokeWidth( vg, 1.0 );
  nvgStroke( vg );
}

void draw_donut( NVGcontext *vg, const vec2_t center, const float r1, const float r2 )
{
  nvgBeginPath( vg );
  nvgCircle( vg, center.x, center.y, r1 );
  nvgCircle( vg, center.x, center.y, r2 );
  nvgStrokeColor( vg, nvgRGBA(64, 64, 64, 255));
  nvgStrokeWidth( vg, 1 );
  nvgStroke( vg );
}

void draw_points( NVGcontext *vg, const vec2_t* pts, const size_t n_pts )
{
  float size = 1.0f;
  nvgBeginPath( vg );
  for( size_t i = 0; i < n_pts; i++ )
  {
    nvgCircle( vg, pts[i].x, pts[i].y, size );
  }
  nvgFillColor( vg, nvgRGBA( 0, 0, 0, 255 ) );
  nvgFill( vg );
}

void draw_point( NVGcontext *vg, const vec2_t pt )
{
  float size = 1.0f;
  nvgBeginPath( vg );
  nvgCircle( vg, pt.x, pt.y, size );
  nvgFillColor( vg, nvgRGBA( 255, 0, 0, 255 ) );
  nvgFill( vg );
}

//=================================================================================================

int main()
{
  // Sampling params
  float radius = 10.0f;
  int k = 10;
  vec2_t domain_size = {.x = 512.0f, .y = 512.0f };

  float cell_size = radius / sqrt(2);
  domain_size.x = ((int)domain_size.x / (int)cell_size) * cell_size;
  domain_size.y = ((int)domain_size.y / (int)cell_size) * cell_size;

  srand( 12345ULL );
  GLFWwindow* window;
  NVGcontext* vg = NULL;
  double prevt = 0;

  if (!glfwInit()) {
    printf("Failed to init GLFW.");
    return -1;
  }

  glfwWindowHint(GLFW_OPENGL_DEBUG_CONTEXT, 1);
  window = glfwCreateWindow( domain_size.x + 2*PADDING, domain_size.y + 2*PADDING, "Poisson Disk Sampling", NULL, NULL);

  if (!window) {
    glfwTerminate();
    return -1;
  }

  glfwMakeContextCurrent(window);

  glewExperimental = GL_TRUE;
  if(glewInit() != GLEW_OK) {
    printf("Could not init glew.\n");
    return -1;
  }
  // GLEW generates GL error because it calls glGetString(GL_EXTENSIONS), we'll consume it here.
  glGetError();

  vg = nvgCreateGL3(NVG_ANTIALIAS | NVG_STENCIL_STROKES | NVG_DEBUG);

  if (vg == NULL) {
    printf("Could not init nanovg.\n");
    return -1;
  }

  glfwSwapInterval(0);
  glfwSetTime(0);
  prevt = glfwGetTime();

  // vec2_t center = {.x = 156, .y = 156 };
  // vec2_t *pts = generate_circular_samples_range( center, radius, 2*radius, k );

  // int n_pts = ((int)domain_size.x / (int)cell_size) * ((int)domain_size.y / (int)cell_size);
  // vec2_t *pts = random_sampling( n_pts, domain_size );

  // int n_pts = 0;
  // vec2_t *pts = stratified_sampling( domain_size, cell_size, &n_pts );

  int n_pts = 0;
  vec2_t *pts = poisson_disk_sampling( domain_size, radius, k, &n_pts );
  test_poisson_samples( pts, n_pts, radius );

  while (!glfwWindowShouldClose(window))
  {
    double t = glfwGetTime();
    double dt = t - prevt;
    prevt = t;

    int winWidth, winHeight;
    int fbWidth, fbHeight;
    float pxRatio;

    // Calculate pixel ration for hi-dpi devices.
    glfwGetWindowSize(window, &winWidth, &winHeight);
    glfwGetFramebufferSize(window, &fbWidth, &fbHeight);

    pxRatio = (float)fbWidth / (float)winWidth;

    // Update and render
    glViewport(0, 0, fbWidth, fbHeight);
    glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT|GL_STENCIL_BUFFER_BIT);

    nvgBeginFrame(vg, winWidth, winHeight, pxRatio);
    
    // draw_grid( vg, domain_size, cell_size );

    draw_points( vg, pts, n_pts );

    // draw_point( vg, center );
    // draw_donut( vg, center, radius, 2*radius );
    // draw_points( vg, pts, k );

    nvgEndFrame(vg);

    glfwSwapBuffers(window);
    glfwPollEvents();
  }

  nvgDeleteGL3(vg);
  glfwTerminate();
  return 0;
}