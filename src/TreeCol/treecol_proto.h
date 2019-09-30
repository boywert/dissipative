#ifdef TREE_RAD
#define PROJECT_COLUMN project_column_
#define GET_ANGULAR_COORDS get_angular_coords_
#define GET_PIXELS_FOR_XYZ_AXES get_pixels_for_xyz_axes_
#define BUILDCOL2D  buildcol2d_
#define CALCULATE_PIXEL_CENTRES calculate_pixel_centres_
#define CREATETRIGLOOKUP createtriglookup_


void PROJECT_COLUMN(double* column, double* columnH2, double* columnCO, double* dx_heal, double* dy_heal, double* dz_heal, double* ang_radius, 
                    double projection[NPIX], double projectionH2[NPIX], double projectionCO[NPIX],int* id);

void GET_ANGULAR_COORDS(int* k, double* theta, double* phi);
void GET_PIXELS_FOR_XYZ_AXES(int pixels[6]);
void BUILDCOL2D(void);
void CALCULATE_PIXEL_CENTRES(void);
void CREATETRIGLOOKUP(void);
#endif
