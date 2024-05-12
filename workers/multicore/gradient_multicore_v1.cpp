#include <stdio.h> 
#include <stdlib.h> 
#include <dirent.h> 
#include <signal.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "mkl.h"

#include <zmq.h>
#include <czmq.h>
#include <unistd.h>
#include "cJSON.h"

static volatile int running = 1;

void interrupt_handler(int signal) {
    running = 0;
}

double get_time_ms() {
    struct timespec ts;
    if (clock_gettime(CLOCK_MONOTONIC, &ts) == 0) {
        return ts.tv_sec * 1000 + ts.tv_nsec / 1000000.0;
    } else {
        return 0;
    }
}

double get_time_s() {
    struct timespec ts;
    if (clock_gettime(CLOCK_MONOTONIC, &ts) == 0) {
        return ts.tv_sec + ts.tv_nsec / 1e9;
    } else {
        return 0;
    }
}

int compare_strings(const void *va, const void *vb) {

    char **a = (char **) va;
    char **b = (char **) vb;
    return strcmp(*a, *b);
}

int read_files(const char *path, const char *pattern, char **files) {

    int n_files = 0;
    DIR *d = opendir(path);
    struct dirent *dir;
    if (d) {
        while ((dir = readdir(d)) != NULL) {
            if (strstr(dir->d_name, pattern) != NULL) {
                files[n_files] = (char *) malloc(1000 * sizeof(char));
                strcpy(files[n_files], path);
                if (path[strlen(path) - 1] != '/') {
                    strcat(files[n_files], (char *) "/");
                }
                strcat(files[n_files], dir->d_name);
                n_files++;
            }
        }
        closedir(d);
    }

    qsort(files, n_files, sizeof(char *), compare_strings);
    return n_files;
}

struct SparseMatrix {
    int n_nz;
    int n_rows;
    int n_cols;
    int *rows;
    int *cols;
    double *vals;
    int *rows_t;
    int *cols_t;
    double *vals_t;

    void malloc_cpu() {
        rows = (int *) malloc(n_nz*sizeof(int));
        cols = (int *) malloc(n_nz*sizeof(int));
        vals = (double *) malloc(n_nz*sizeof(double));
        n_rows = 0;
        n_cols = 0;
    }

    void free_cpu() {
        free(rows);
        free(cols);
        free(vals);
        free(rows_t);
        free(cols_t);
        free(vals_t);
    }
};

struct Region {
    char *name;
    int n_voxels;
    double min;
    double max;
    double avg;

    double f; // Objective function evaluation
    double eud;
    double dF_dEUD;
    double sum_alpha;
    // Virtual EUD to control PTV overdosage
    // Hardcoded to eud + 1 for now
    double v_f;
    double v_eud;
    double v_dF_dEUD;
    double v_sum_alpha;

    bool is_optimized;
    bool is_ptv;
    double pr_min;
    double pr_max;
    double pr_avg_min;
    double pr_avg_max;
    double *grad_avg;
    
    double pr_eud;
    int penalty;
    int alpha;

    void set_targets(bool t_ptv, double t_min, double t_avg_min, double t_avg_max, double t_max, 
                     double t_eud, int t_alpha, int t_penalty) {
        if (t_eud < 0 && t_min < 0 && t_max < 0 && 
            t_avg_min < 0 && t_avg_max < 0) {
            is_optimized = false;
        } else {
            is_optimized = true;
            is_ptv = t_ptv;
            pr_min = t_min;
            pr_max = t_max;
            pr_avg_min = t_avg_min;
            pr_avg_max = t_avg_max;
            pr_eud = t_eud;
            alpha = t_alpha;
            penalty = t_penalty;
            f = 0;
            v_f = 0;
            eud = 0;
            v_eud = 0;
            dF_dEUD = 0;
            v_dF_dEUD = 0;
            sum_alpha = 0;
            v_sum_alpha = 0;
        }
    }
};

#define BEAM_MAP_X 120
#define BEAM_MAP_Y 120

struct Plan {
    char *name;
    int n_beams;
    int n_beamlets;
    int *n_beamlets_beam;
    int n_voxels;
    int n_regions;
    int n_plans;
    double dose_grid_scaling;
    Region* regions;
    char *voxel_regions;
    SparseMatrix spm;
    double *fluence;
    double *smoothed_fluence;
    double *doses;
    char *files[100];
    int *beam_maps;
    sparse_matrix_t m;
    sparse_matrix_t m_t;
    struct matrix_descr descr;

    double *f0;
    double *f1;

    void check_line(int result) {
        if (result < 0) {
            fprintf(stderr, "ERROR in %s (%s:%d): Unable to read line.\n", 
                    __func__, __FILE__, __LINE__);
            exit(EXIT_FAILURE);
        }
    }

    char* get_str(char *line, size_t len) {
        char *temp = (char *) malloc(len*sizeof(char));
        snprintf(temp, len, "%s", line);
        temp[strcspn(temp, "\r\n")] = 0; // Remove newline
        return temp;
    }

    int get_int(char *line, char **end) {
        return strtoll(line, end, 10);
    }

    float get_float(char *line, char **end) {
        return strtof(line, end);
    }

    void parse_config(const char *path) {
        int n_files = read_files(path, "m_", files);

        FILE *f = fopen(files[0], "r");
        if (f == NULL) {
            fprintf(stderr, "ERROR in %s (%s:%d): Unable to open file.\n", 
                    __func__, __FILE__, __LINE__);
            exit(EXIT_FAILURE);
        }
        printf("Loading %s...\n", files[0]);

        char *line = NULL;
        char *end = NULL;
        size_t len = 0;

        check_line(getline(&line, &len, f));
        name = get_str(line, len);

        check_line(getline(&line, &len, f));
        n_beams = get_int(line, &end);

        n_beamlets = 0;
        n_beamlets_beam = (int *)malloc(n_beams * sizeof(int));
        for (int i = 0; i < n_beams; i++) {
            check_line(getline(&line, &len, f));
            int index = get_int(line, &end);
            int beamlets = get_int(end, &line);
            n_beamlets_beam[index - 1] = beamlets;
            n_beamlets += beamlets;
        }

        check_line(getline(&line, &len, f));
        n_voxels = get_int(line, &end);

        check_line(getline(&line, &len, f));
        dose_grid_scaling = get_float(line, &end);

        check_line(getline(&line, &len, f));
        n_regions = get_int(line, &end);

        regions = (Region *) malloc(n_plans*n_regions*sizeof(Region));
        for (int i = 0; i < n_regions; i++) {
            check_line(getline(&line, &len, f));
            get_int(line, &end);
            char *name = get_str(end + 1, len);
            regions[i].name = name;
            regions[i].n_voxels = 0;
        }
        for (int k = 1; k < n_plans; k++) {
            for (int i = 0; i < n_regions; i++) {
                regions[k*n_regions + i].name = regions[i].name;
                regions[k*n_regions + i].n_voxels = 0;
            }
        }

        line = NULL;
        len = 0;
        while (getline(&line, &len, f) != -1) {
            fprintf(stderr, "[WARNING] Line not processed: %s", line);
        }

        fclose(f);
        free(files[0]);
    }

    void parse_voxel_regions(const char *path) {
        int n_files = read_files(path, "v_", files);

        FILE *f = fopen(files[0], "r");
        if (f == NULL) {
            fprintf(stderr, "ERROR in %s (%s:%d): Unable to open file.\n",
                    __func__, __FILE__, __LINE__);
            exit(EXIT_FAILURE);
        }
        printf("Loading %s...\n", files[0]);

        voxel_regions = (char *) malloc(n_voxels*n_regions*sizeof(char));
        char line[1024];
        int num = 0;
        int offset = 0;
        while (fgets(line, sizeof line, f)) {
            if (sscanf(line, "%d", &num)) {
                for (int i = 0; i < n_regions; i++) {
                    voxel_regions[offset + i*n_voxels] = num & 1;
                    num >>= 1;
                }
                offset++;
            } else {
                fprintf(stderr, "ERROR in %s (%s:%d): Unable to read voxel regions.\n",
                        __func__, __FILE__, __LINE__);
                exit(EXIT_FAILURE);
            }
        }

        for (int i = 0; i < n_regions; i++) {
            for (int j = 0; j < n_voxels; j++) {
                if (voxel_regions[i*n_voxels + j]) {
                    for (int k = 0; k < n_plans; k++) {
                        regions[k*n_regions + i].n_voxels += 1;
                    }
                }
            }
        }

        fclose(f);
        free(files[0]);
    }

    void load_spm(const char *path) {
        int n_files = read_files(path, "d_", files);

        FILE **fp = (FILE **) malloc(n_files*sizeof(FILE *));
        for (int i = 0; i < n_files; i++) {
            fp[i] = fopen(files[i], "r");
            if (fp[i] == NULL) {
                fprintf(stderr, "ERROR in %s (%s:%d): Unable to open file.\n",
                        __func__, __FILE__, __LINE__);
                exit(EXIT_FAILURE);
            }
            int n_nz = 0;
            int count = fscanf(fp[i], "%d", &n_nz);
            spm.n_nz += n_nz;
        }

        spm.malloc_cpu();

        int idx = 0;
        int offset = 0;
        for (int i = 0; i < n_files; i++) {
            printf("Loading %s... ", files[i]);

            int n_read = 0;
            while (true) {
                int row, col;
                double val;
                int count = fscanf(fp[i], "%d %d %lf", &row, &col, &val);
                if(count == EOF || !count) {
                    break;
                }

                int new_col = offset + col;
                spm.rows[idx] = row;
                spm.cols[idx] = new_col;
                spm.vals[idx] = val;
                idx++;
                n_read++;

                if (row > spm.n_rows) {
                    spm.n_rows = row;
                }
                if (new_col > spm.n_cols) {
                    spm.n_cols = new_col;
                }
            }

            printf("%d values read.\n", n_read);
            offset = spm.n_cols + 1;
            fclose(fp[i]);
            free(files[i]);
        }

        spm.n_rows++;
        // Sometimes there's missing voxels,
        // but we want the dimensions to match for SpMM
        if (spm.n_rows < n_voxels) {
            spm.n_rows = n_voxels;
        }
        spm.n_cols++;

        free(fp);

        descr.type = SPARSE_MATRIX_TYPE_GENERAL;
        mkl_sparse_d_create_coo(&m, SPARSE_INDEX_BASE_ZERO, spm.n_rows, spm.n_cols, spm.n_nz, spm.rows, spm.cols, spm.vals);
        mkl_sparse_convert_csr(m, SPARSE_OPERATION_TRANSPOSE, &m_t);
        mkl_sparse_convert_csr(m, SPARSE_OPERATION_NON_TRANSPOSE, &m);
        mkl_sparse_set_mv_hint(m, SPARSE_OPERATION_NON_TRANSPOSE, descr, 1e6);
        mkl_sparse_set_mv_hint(m_t, SPARSE_OPERATION_NON_TRANSPOSE, descr, 1e6);
        mkl_sparse_optimize(m);
        mkl_sparse_optimize(m_t);
    
    }

    void load_fluence(const char *path, const char *prefix) {
        int n_files = read_files(path, prefix, files);

        FILE **fp = (FILE **) malloc(n_files*sizeof(FILE *));
        int idx = 0;
        for (int i = 0; i < n_files; i++) {
            printf("Loading %s... ", files[i]);
            fp[i] = fopen(files[i], "r");
            if (fp[i] == NULL) {
                fprintf(stderr, "ERROR in %s (%s:%d): Unable to open file.\n",
                        __func__, __FILE__, __LINE__);
                exit(EXIT_FAILURE);
            }

            int n_read = 0;
            while (true) {
                int count = fscanf(fp[i], "%lf", &(fluence[idx]));
                if(count == EOF || !count) {
                    break;
                }

                idx++;
                n_read++;
            }

            printf("%d values read.\n", n_read);
            fclose(fp[i]);
            free(files[i]);
        }

        free(fp);

        for (int i = 1; i < n_plans; i++) {
            for (int j = 0; j < n_beamlets; j++) {
                fluence[i*n_beamlets + j] = fluence[j];
            }
        }

    }

    void load_coords(const char *path) {
        int n_files = read_files(path, "xcoords_", files);

        beam_maps = (int *) malloc(n_beams*BEAM_MAP_Y*BEAM_MAP_X*sizeof(int));

        for (int i = 0; i < n_beams*BEAM_MAP_Y*BEAM_MAP_X; i++) {
            beam_maps[i] = -1;
        }

        int idx = 0;
        FILE **fp = (FILE **) malloc(n_files*sizeof(FILE *));
        for (int i = 0; i < n_files; i++) {
            fp[i] = fopen(files[i], "r");
            if (fp[i] == NULL) {
                fprintf(stderr, "ERROR in %s (%s:%d): Unable to open file.\n",
                        __func__, __FILE__, __LINE__);
                exit(EXIT_FAILURE);
            }
            printf("Loading %s... ", files[i]);

            char ignored[1024];
            int n_read = 0;
            while (true) {
                int col, row;
                int count = fscanf(fp[i], "%d %d", &col, &row);
                if(count == EOF) {
                    break;
                } else if (!count) {
                    fgets(ignored, sizeof(ignored), fp[i]);
                } else if (count == 1) {
                    // Header values, ignored
                    continue;
                } else if (count == 2) {
                    beam_maps[i*BEAM_MAP_Y*BEAM_MAP_X + col*BEAM_MAP_Y + row] = idx;
                    n_read++;
                    idx++;
                } else {
                    fprintf(stderr, "ERROR in %s (%s:%d): While reading coordinate file.\n",
                            __func__, __FILE__, __LINE__);
                    exit(EXIT_FAILURE);
                }
            }
            printf("%d coordinates read.\n", n_read);
            fclose(fp[i]);
            free(files[i]);
        }


        free(fp);
    }

    void init_fluence(float value) {
        for (int i = 0; i < n_plans; i++) {
            for (int j = 0; j < n_beamlets; j++) {
                fluence[i*n_beamlets + j] = value;
            }
        }
    }

    void print() {
        printf("Name: %s\n", name);
        printf("Number of beams: %d\n", n_beams);
        for (int i = 0; i < n_beams; i++) {
            printf("  Beam %d: %d beamlets\n", i + 1, n_beamlets_beam[i]);
        }
        printf("Total: %d beamlets\n", n_beamlets);
        printf("Number of voxels: %d\n", n_voxels);
        printf("Dose Grid Scaling: %e\n", dose_grid_scaling);
        printf("Number of regions: %d\n", n_regions);
        for (int i = 0; i < n_regions; i++) {
            printf("  Region %2d: %-16s %8d voxels\n", i, regions[i].name, regions[i].n_voxels);
        }
        printf("Dose matrix: %d x %d with %d nonzeros.\n", spm.n_rows, spm.n_cols, spm.n_nz);
        printf("Number of plans: %d.\n", n_plans);
    }

    void compute_dose() {
        memset(doses, 0, n_plans*n_voxels*sizeof(*doses));

        double alpha = 1.0, beta = 0.0;

        double start_time = get_time_s();
        mkl_sparse_d_mm(SPARSE_OPERATION_NON_TRANSPOSE, alpha, m, descr, SPARSE_LAYOUT_COLUMN_MAJOR, fluence, n_plans, n_beamlets, beta, doses, n_voxels);
        double elapsed = get_time_s() - start_time;
        //printf("Compute dose mm: %.4f seconds.\n", elapsed);

        #pragma omp parallel for
        for (int i = 0; i < n_voxels*n_plans; i++) {
            doses[i] *= dose_grid_scaling;
        }
    }

    void stats() {
        #pragma omp parallel for collapse(2)
        for (int k = 0; k < n_plans; k++) {
            for (int i = 0; i < n_regions; i++) {
                Region *r = &regions[k*n_regions + i];

                r->min = 1e10;
                r->max = 0;
                r->avg = 0;
                r->sum_alpha = 0;
                r->v_sum_alpha = 0;
                for (int j = 0; j < n_voxels; j++) {
                    if (voxel_regions[i*n_voxels + j]) {
                        double dose = doses[k*n_voxels + j];
                        if (r->min > dose) {
                            r->min = dose;
                        }
                        if (r->max < dose) {
                            r->max = dose;
                        }
                        r->avg += dose;

                        if (dose > 0) {
                            r->sum_alpha += pow((double) dose, 
                                                (double) r->alpha);
                            if (r->is_ptv) {
                                r->v_sum_alpha += pow((double) dose, 
                                                      (double) -r->alpha);
                            }
                        }
                    }
                }
                r->avg /= r->n_voxels;
                r->eud = pow(r->sum_alpha/r->n_voxels, 
                             1.0/r->alpha);
                if (r->is_ptv) {
                    r->v_eud = pow(r->v_sum_alpha/r->n_voxels, 
                                   1.0/-r->alpha);
                }
                if (r->is_optimized) {
                    int n = r->penalty;
                    int pd = r->pr_eud;
                    double eud = r->eud;
                    int v_pd = pd + 1; // Hardcoded virtual PTV prescribed dose
                    double v_eud = r->v_eud;
                    if (r->is_ptv) {
                        r->f = 1/(1 + pow(pd/eud, n));
                        r->dF_dEUD =  (n*r->f/eud) * pow(pd/eud, n);
                        // Virtual EUD to control PTV over-dosage
                        r->v_f = 1/(1 + pow(v_eud/v_pd, n));
                        r->v_dF_dEUD = -(n*r->v_f/v_eud) * pow(v_eud/v_pd, n);
                    } else {
                        r->f = 1/(1 + pow(eud/pd, n));
                        r->dF_dEUD = -(n*r->f/eud) * pow(eud/pd, n);
                    }
                }
            }
        }
    }

    void print_table(int pid) {
        printf("%2d    Region          Min         Avg         Max         EUD     dF_dEUD         v_EUD     v_dF_dEUD      f       v_f\n", pid); 
        for (int i = 0; i < n_regions; i++) {
            Region r = regions[pid*n_regions + i];
            if (r.is_optimized) {
                printf("%-15s %11.6lf %11.6lf %11.6lf %11.6lf %11.6lf %11.6lf %11.6lf %11.6lf %11.6lf\n", r.name, r.min, r.avg, r.max, r.eud, r.dF_dEUD, r.v_eud, r.v_dF_dEUD, r.f, r.v_f);
            }
        }
    }

    void load(const char *plan_path, const char *fluence_path, const char *fluence_prefix) {
        parse_config(plan_path);
        parse_voxel_regions(plan_path);
        load_spm(plan_path);

        fluence = (double *) malloc(n_plans*n_beamlets*sizeof(double));
        smoothed_fluence = (double *) malloc(n_plans*n_beamlets*sizeof(double));
        doses = (double *) malloc(n_plans*n_voxels*sizeof(double));

        f0 = (double *) malloc(n_plans*sizeof(double));
        f1 = (double *) malloc(n_plans*sizeof(double));

        load_coords(plan_path);

        load_fluence(fluence_path, fluence_prefix);
        //init_fluence(1e-2);
        print();
    }


    void smooth_cpu() {
        int n_neighbors = 8;
        int sum_weights = 1000;

        #pragma omp parallel for collapse(2)
        for (int k = 0; k < n_plans; k++) {
            for (int i = 0; i < n_beams; i++) {
                int *neighbors = (int *) malloc(n_neighbors*sizeof(int));
                int poff = k*n_beamlets;
                for (int y = 0; y < BEAM_MAP_Y; y++) {
                    for (int x = 0; x < BEAM_MAP_X; x++) {
                        int offset = i*BEAM_MAP_Y*BEAM_MAP_X;
                        int idx = beam_maps[offset + BEAM_MAP_X*y + x];
                        float center_weight = sum_weights - n_neighbors;
                        if (idx >= 0 && idx < n_beamlets) {
                            smoothed_fluence[poff + idx] = 0;
                            neighbors[0] = beam_maps[offset + BEAM_MAP_Y*(y-1) + (x-1)];
                            neighbors[1] = beam_maps[offset + BEAM_MAP_Y*(y-1) + (x  )];
                            neighbors[2] = beam_maps[offset + BEAM_MAP_Y*(y-1) + (x+1)];
                            neighbors[3] = beam_maps[offset + BEAM_MAP_Y*(y  ) + (x-1)];
                            neighbors[4] = beam_maps[offset + BEAM_MAP_Y*(y  ) + (x+1)];
                            neighbors[5] = beam_maps[offset + BEAM_MAP_Y*(y+1) + (x-1)];
                            neighbors[6] = beam_maps[offset + BEAM_MAP_Y*(y+1) + (x  )];
                            neighbors[7] = beam_maps[offset + BEAM_MAP_Y*(y+1) + (x+1)];

                            //if (neighbors[3] < 0 || neighbors[4] < 0) {
                            //    // This is a border beamlet, ignore other rows
                            //    neighbors[0] = -1;
                            //    neighbors[1] = -1;
                            //    neighbors[2] = -1;
                            //    neighbors[5] = -1;
                            //    neighbors[6] = -1;
                            //    neighbors[7] = -1;

                            //}
                            for (int j = 0; j < n_neighbors; j++) {
                                if (neighbors[j] >= 0 && neighbors[j] < n_beamlets) {
                                    smoothed_fluence[poff + idx] += fluence[poff + neighbors[j]];
                                } else {
                                    center_weight += 0.8;
                                }
                            }
                            smoothed_fluence[poff + idx] += center_weight*fluence[poff + idx];
                            smoothed_fluence[poff + idx] /= sum_weights;
                        }
                    }
                }
                free(neighbors);
            }
        }

        #pragma omp parallel for
        for (int i = 0; i < n_beamlets*n_plans; i++) {
            fluence[i] = smoothed_fluence[i];
        }
    }
};

void voxels_eud(Plan plan, int rid, int pid, double *voxels) {
    Region *r = &plan.regions[pid*plan.n_regions + rid];

    #pragma omp parallel for
    for (int i = 0; i < plan.n_voxels; i++) {
        double dose = plan.doses[pid*plan.n_voxels + i];
        if (plan.voxel_regions[rid*plan.n_voxels + i]) {
            double dEUD_dd = r->eud*pow(dose, r->alpha - 1)/r->sum_alpha;
            voxels[i] = r->dF_dEUD * dEUD_dd;
            if (r->is_ptv) {
                dEUD_dd = r->v_eud*pow(dose, -r->alpha - 1)/r->v_sum_alpha;
                voxels[plan.n_voxels + i] = r->v_dF_dEUD * dEUD_dd;
            }
        } else {
            voxels[i] = 0;
            if (r->is_ptv) {
                voxels[plan.n_voxels + i] = 0;
            }
        }
    }
}

double penalty(Plan plan, unsigned int pid) {
    double penalty = 0;

    for (int i = 0; i < plan.n_regions; i++) {
        Region region = plan.regions[pid*plan.n_regions + i];
        if (region.is_optimized) {
            if (region.pr_min > 0 &&
                region.min < region.pr_min) {
                penalty += region.pr_min - region.min;
            }
            if (region.pr_max > 0 && 
                region.max > region.pr_max) {
                penalty += region.max - region.pr_max;
            }
            if (region.pr_avg_min > 0 && 
                region.avg < region.pr_avg_min) {
                penalty += region.pr_avg_min - region.avg;
            }
            if (region.pr_avg_max > 0 && 
                region.avg > region.pr_avg_max) {
                penalty += region.avg - region.pr_avg_max;
            }
        }
    }
    return penalty;
}

double objective(Plan plan, unsigned int pid) {
    double objective = 1;

    for (int i = 0; i < plan.n_regions; i++) {
        Region region = plan.regions[pid*plan.n_regions + i];
        if (region.is_optimized) {
            objective *= region.f;
            if (region.is_ptv) {
                objective *= region.v_f;
            }
        }
    }
    return objective;
}

void vector_stats(const char *name, double *vector, int n_values) {
    double min = 1e10, max = 0, avg = 0;
    for (int i = 0; i < n_values; i++) {
        if (vector[i] < min) {
            min = vector[i];
        }
        if (vector[i] > max) {
            max = vector[i];
        }
        avg += vector[i];
    }
    avg /= n_values;

    printf("%s: %f %f %f\n", name, min, max, avg);
}

void reduce_gradient(double *voxels, int n_voxels, int n_gradients, int n_plans) {
    for (int k = 0; k < n_plans; k++) {
        unsigned int poff = k*n_voxels*n_gradients;
        for (int i = 0; i < n_voxels; i++) {
            double influence = 0;
            for (int j = 0; j < n_gradients; j++) {
                influence += voxels[poff + j*n_voxels + i];
            }
            voxels[k*n_voxels + i] = influence;
        }
    }
}

void apply_gradient(double *gradient, double *momentum, int n_beamlets, float step, double *fluence, int n_plans) {

    int beta = 0.9;

    #pragma omp parallel for
    for (int i = 0; i < n_plans*n_beamlets; i++) {
        momentum[i] = beta*momentum[i] + (1-beta)*gradient[i];
        fluence[i] += step*momentum[i];

        if (fluence[i] < 0) {
            fluence[i] = 0;
        }
        if (fluence[i] > 0.3) {
            fluence[i] = 0.3;
        }
    }
}

int descend(Plan plan, double *voxels, double *gradient, double *momentum, float step) {
    memset(voxels, 0, plan.n_plans*plan.n_voxels*sizeof(*voxels));
    memset(gradient, 0, plan.n_plans*plan.n_beamlets*sizeof(*gradient));

    int n_gradients = 0;
    int offset = 0;

    for (int k = 0; k < plan.n_plans; k++) {
        for (int i = 0; i < plan.n_regions; i++) {
            Region *region = &plan.regions[i];
            if (region->is_optimized) {
                voxels_eud(plan, i, k, &voxels[offset]);
                offset += plan.n_voxels;
                n_gradients++;
                if (region->is_ptv) {
                    offset += plan.n_voxels;
                    n_gradients++;
                }
            }
        }
    }

    n_gradients /= plan.n_plans;
    reduce_gradient(voxels, plan.n_voxels, n_gradients, plan.n_plans);
    
    double alpha = 1.0, beta = 0.0;
    double start_time = get_time_s();
    //mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, alpha, plan.m_t, plan.descr, voxels, beta, gradient);
    mkl_sparse_d_mm(SPARSE_OPERATION_NON_TRANSPOSE, alpha, plan.m_t, plan.descr,SPARSE_LAYOUT_COLUMN_MAJOR, 
                    voxels, plan.n_plans, plan.n_voxels, beta, gradient, plan.n_beamlets);
    double elapsed = get_time_s() - start_time;
    //printf("Descend mm: %.4f seconds.\n", elapsed);
    
    apply_gradient(gradient, momentum, plan.n_beamlets, step, plan.fluence, plan.n_plans);

    return n_gradients;
}

void optimize(Plan plan) {

    int gradients_per_region = 3; // Warning, hardcoded!
    // *2 because we need two for PTVs, but we're wasting space on organs.
    double *voxels = (double *) malloc(plan.n_plans*plan.n_voxels*plan.n_regions*2*sizeof(*voxels)); 
    double *gradient = (double *) malloc(plan.n_plans*plan.n_beamlets*sizeof(*gradient));
    double *momentum = (double *) malloc(plan.n_plans*gradients_per_region*plan.n_regions*plan.n_beamlets*sizeof(*momentum));
    memset(momentum, 0, plan.n_plans*gradients_per_region*plan.n_regions*plan.n_beamlets*sizeof(*momentum));

    plan.compute_dose();
    plan.stats();
    //printf("Initial solution:\n");
    //for (int k = 0; k < plan.n_plans; k++) {
    //    plan.print_table(k);
    //}

    int rid_sll = 3;
    int rid_slr = 4;

    double step = 2e-7;
    double decay = 1e-7;
    double min_step = 1e-1;
    double start_time = get_time_s();
    double current_time;

    int it = 0;
    while (running && get_time_s() - start_time < 600000) {
        descend(plan, voxels, gradient, momentum, step);
        plan.smooth_cpu();
        plan.compute_dose();
        plan.stats();
        //break;

        if (it % 100 == 0) {
            current_time = get_time_s();

            printf("\n[%.3f] Iteration %d %e\n", current_time - start_time, it, step);

            for (int k = 0; k < plan.n_plans; k++) {
                unsigned int poff = k*plan.n_regions;

                double pen = penalty(plan, k);
                double obj = plan.regions[poff + rid_sll].avg + plan.regions[poff + rid_slr].avg;
                double obj2 = objective(plan, k);
                printf("%2d penalty: %9.6f\n", k, pen);
                printf("%2d    obj: %9.6f\n", k, obj); 
                printf("%2d   obj2: %9.24f\n", k, obj2);
                plan.print_table(k);
            }
        }
        if (step > min_step) 
            step = step/(1 + decay*it);
        it++;
        if (it == 2000)
            break;
    }

    double elapsed = get_time_s() - start_time;
    printf("\nRan %d iterations in %.4f seconds (%.4f ms/it) \n", it, elapsed, elapsed*1000/it/plan.n_plans);
    for (int k = 0; k < plan.n_plans; k++) {
        unsigned int poff = k*plan.n_regions;

        double pen = penalty(plan, k);
        double obj = plan.regions[poff + rid_sll].avg + plan.regions[poff + rid_slr].avg;
        printf("[%2d] penalty (f0): %9.6f\n", k, pen);
        printf("[%2d]     obj (f1): %9.6f\n", k, obj); 
        plan.f0[k] = pen;
        plan.f1[k] = obj;
        //plan.print_table(k);
    }

    free(voxels);
    free(gradient);
    free(momentum);
}

void json_parse_roi(cJSON *json, const char *name, Region *region) {
    cJSON *roi = cJSON_GetObjectItem(json, name);
    if (roi) {
        printf("%16s: %s\n",name, cJSON_PrintUnformatted(roi));
        cJSON *roi_a = cJSON_GetObjectItem(roi, "a");
        cJSON *roi_eud = cJSON_GetObjectItem(roi, "EUD_0");
        cJSON *roi_n = cJSON_GetObjectItem(roi, "n");
        if (roi_a) {
            printf("%s: %s\n", "a", cJSON_PrintUnformatted(roi_a));
            region->alpha = roi_a->valuedouble;
            //cJSON_Delete(roi_a);
        }
        if (roi_n) {
            printf("%s: %s\n", "n", cJSON_PrintUnformatted(roi_n));
            region->penalty = roi_n->valuedouble;
            //cJSON_Delete(roi_n);
        }
        if (roi_eud) {
            printf("%s: %s\n", "EUD_0", cJSON_PrintUnformatted(roi_eud));
            region->pr_eud = roi_eud->valuedouble;
            //cJSON_Delete(roi_eud);
        }
        

         //cJSON_Delete(roi);
    }
   
    
}

int main(int argc, char **argv) {

    signal(SIGINT, interrupt_handler);

    int plan_n = atoi(argv[1]);
    const char* plan_path = argv[2];
    const char* out_path = argv[3];
    const char* fluence_path;
    const char* fluence_prefix;

    if (argc > 4) {
        fluence_path = argv[4];
        fluence_prefix = argv[5];
    } else {
        // We use the starting plan from Eclipse
        fluence_path = plan_path;
        fluence_prefix = "x_PARETO";
    }


    //int n_plans = atoi(argv[6]);


	zsock_t *worker = zsock_new(ZMQ_REQ);
    char hostname[4];
    char identity[10];
    gethostname(hostname, 3);
    sprintf(identity, "%s-%04d", hostname, getpid());
    zmq_setsockopt(zsock_resolve(worker), ZMQ_IDENTITY, identity, strlen(identity));
	zsock_connect(worker, "tcp://bullxual:3560");

	zframe_t *frame = zframe_new("READY", 5);
	zframe_send(&frame, worker, 0);
    printf("Worker %s started.\n", identity);


    while (true) {
        zmsg_t *msg = zmsg_recv(worker);
        if (!msg) break; // Error!
        char *data = zframe_strdup(zmsg_last(msg));
        printf("Received raw data: %s\n", data);
        
        cJSON *json = cJSON_Parse(data);
        cJSON *json_plan = NULL;
        cJSON *json_plans = cJSON_GetObjectItemCaseSensitive(json, "Plans");

        int n_plans = 0;
        cJSON_ArrayForEach(json_plan, json_plans) {
            n_plans++;
        }
        printf("Received number of plans: %d\n", n_plans);

        // Esto debería hacerse antes de recibir el trabajo, quitar la dependencia de n_plans en el load.
        Plan plan = {};
        plan.n_plans = n_plans;
        plan.load(plan_path, fluence_path, fluence_prefix);

        if (plan_n == 5) {
            for (int i = 0; i < plan.n_plans; i++) {
                plan.regions[i*plan.n_regions +  0].set_targets(false,    -1,    -1,    -1,    -1,    -1,  10,   5);
                plan.regions[i*plan.n_regions +  1].set_targets(false,    -1,    -1,    -1, 74.25, 74.25,  40,   5);
                plan.regions[i*plan.n_regions +  2].set_targets(false,    -1,    -1,    -1,    70,    70,  10,   5);
                plan.regions[i*plan.n_regions +  3].set_targets(false,    -1,    -1,    26,    -1,  4.76, 1.01, 18.25);
                plan.regions[i*plan.n_regions +  4].set_targets(false,    -1,    -1,    26,    -1,  3.79, 1.31, 11.17);
                plan.regions[i*plan.n_regions +  5].set_targets(false,    -1,    -1,    -1,    50,  1.80, 1.33, 12.79);
                plan.regions[i*plan.n_regions +  6].set_targets(false,    -1,    -1,    -1,    -1,    -1,  10,   5);
                plan.regions[i*plan.n_regions +  7].set_targets(false,    -1,    -1,    -1,    60,    60,  10,   5);
                plan.regions[i*plan.n_regions +  8].set_targets( true, 48.60, 52.92, 55.08, 59.40, 54.00, -65.55, 54.11);
                plan.regions[i*plan.n_regions +  9].set_targets( true, 54.00, 58.80, 61.20, 66.00, 60.00, -87.12, 57.66);
                plan.regions[i*plan.n_regions + 10].set_targets( true, 59.40, 64.67, 67.32, 72.60, 66.00, -33.27, 18.62);
                plan.regions[i*plan.n_regions + 11].set_targets(false,    -1,    -1,    -1,    -1,    -1,  10,   5);
            }
        } else {
            return 1;
        }
            
        int i = 0;
        cJSON_ArrayForEach(json_plan, json_plans) {
            printf("Received parameters for Plan %d:\n", i);
            cJSON *json_params = cJSON_GetObjectItemCaseSensitive(json_plan, "Variables");
            json_parse_roi(json_params, "Salivary Gland R", &plan.regions[i*plan.n_regions +  3]);
            json_parse_roi(json_params, "Salivary Gland L", &plan.regions[i*plan.n_regions +  4]);
            json_parse_roi(json_params, "Spinal Cord +3mm", &plan.regions[i*plan.n_regions +  5]);
            json_parse_roi(json_params, "PTV 54",           &plan.regions[i*plan.n_regions +  8]);
            json_parse_roi(json_params, "PTV 60",           &plan.regions[i*plan.n_regions +  9]);
            json_parse_roi(json_params, "PTV 66",           &plan.regions[i*plan.n_regions + 10]);
            printf("\n");
            i++;
            cJSON_Delete(json_params);
        }

        optimize(plan);

        cJSON *evaluations = cJSON_CreateArray();
        for (int i = 0; i < n_plans; i++) {
            cJSON *evaluation = cJSON_CreateObject();
            cJSON_AddNumberToObject(evaluation, "f0", plan.f0[i]);
            cJSON_AddNumberToObject(evaluation, "f1", plan.f1[i]);
            cJSON_AddItemToArray(evaluations, evaluation);
        }
        free(plan.f0);
        free(plan.f1);
        free(plan.fluence);
        free(plan.smoothed_fluence);
        free(plan.doses);
        free(plan.beam_maps);
        free(plan.voxel_regions);
        free(plan.n_beamlets_beam);
        free(plan.name);
        plan.spm.free_cpu();
        for(int i = 0; i < plan.n_regions;i++){
            free(plan.regions[i].name);
        }
        
        free(plan.regions);
        mkl_sparse_destroy(plan.m);
         mkl_sparse_destroy(plan.m_t);
       

        free(data);
        //cJSON_Delete(json);
        //cJSON_Delete(json_plan);
        //cJSON_Delete(json_plans);

        char *response = cJSON_PrintUnformatted(evaluations);
        printf("Replying: %s\n", response);

        zframe_reset(zmsg_last(msg), response, strlen(response));
		zmsg_send(&msg, worker);
        
    }
   
    

	if (frame)
        zframe_destroy(&frame);
	zsock_destroy(&worker);

    return 0;

}
