#include <unistd.h>
#include <string.h>
#include <thread>
#include <mutex>
#include "raytracer.h"

#define THREADS 32

std::thread threads[THREADS];
int thread_number[THREADS];
std::mutex mutex;

bool check_sphere_hit(const std::vector<Sphere>& spheres, const Ray& ray, float t_min, float t_max, Hit& hit) {
    Hit closest_hit;
    bool has_hit = false;
    auto closest_hit_distance = t_max;
    Material material;

    for(std::size_t i = 0; i < spheres.size(); i++) {
        const auto& sphere = spheres[i];
        if(sphere_hit(sphere, ray, t_min, closest_hit_distance, closest_hit)) {
            has_hit = true;
            closest_hit_distance = closest_hit.t;
            material = sphere.material;
        }
    }

    if(has_hit) {
        hit = closest_hit;
        hit.material = material;
    }
    
    return has_hit;
}

Vector3 trace_ray(const Ray& ray, const std::vector<Sphere>& spheres, int depth) {
    if (depth <= 0) {
        return Vector3(0, 0, 0);
    }

    Hit hit;
    if(check_sphere_hit(spheres, ray, 0.001f, FLT_MAX, hit)) {
        Ray outgoing_ray;
        Vector3 attenuation;

        if(metal_scater(hit.material, ray, hit, attenuation, outgoing_ray)) {
            auto ray_color = trace_ray(outgoing_ray, spheres, depth - 1);
            return Vector3(ray_color.x * attenuation.x, ray_color.y * attenuation.y, ray_color.z * attenuation.z);
        }

        return Vector3(0, 0, 0);
    }

    Vector3 unit_direction = unit_vector(ray.direction);
    auto t = 0.5 * (unit_direction.y + 1.0);
    return Vector3(1.0, 1.0, 1.0) * (1.0 - t) + Vector3(0.5, 0.7, 1.0) * t;
}

int width = IMAGE_WIDTH;
int height = IMAGE_HEIGHT;
int samples = NUM_SAMPLES;
int depth = SAMPLE_DEPTH;
int* image_data;
const auto aspect_ratio = (float) width / height;
Camera camera(Vector3(0,1,1), Vector3(0,0,-1), Vector3(0,1,0), aspect_ratio, 90, 0.0f, 1.5f);
Checksum checksum(0,0,0);
std::vector<Sphere> spheres;

void kernel(int thread_num) {
    int i_start = thread_num*height/THREADS;
    int i_end = (thread_num+1)*height/THREADS-1;
    Checksum local_checksum(0,0,0);
    for(int y = i_end; y >= i_start; y--) {
        for(int x = 0; x < width; x++) {
            Vector3 pixel_color(0,0,0);
            for(int s = 0; s < samples; s++) {
                auto u = (float) (x + random_float()) / (width - 1);
                auto v = (float) (y + random_float()) / (height - 1);
                auto r = get_camera_ray(camera, u, v);
                pixel_color += trace_ray(r, spheres, depth);
            }
	    auto output_color = compute_color(local_checksum, pixel_color, samples);

            int pos = ((height - 1 - y) * width + x) * 3;
            image_data[pos] = output_color.r;
            image_data[pos + 1] = output_color.g;
            image_data[pos + 2] = output_color.b;
        }
    }
    mutex.lock();
    checksum += local_checksum;
    mutex.unlock();
}


int main(int argc, char **argv) {


    // This option parsing is not very interesting.
    int no_output = 0;
    char file_name[256] = "render.ppm";
    int c;
    while ((c = getopt(argc, argv, "d:s:r:n:f:")) != -1)
    {
        switch (c)
        {
            case 'd':
                if (sscanf(optarg, "%d", &depth) != 1)
                    goto error;
                break;
            case 's':
                if (sscanf(optarg, "%d", &samples) != 1)
                    goto error;
                break;
            case 'r':
                if (sscanf(optarg, "%dx%d", &width, &height) != 2)
                    goto error;
                break;
            case 'n':
                if (sscanf(optarg, "%d", &no_output) != 1)
                    goto error;
                break;
            case 'f':
                strncpy(file_name, optarg, sizeof(file_name));
                file_name[255] = '\0'; // safe-guard null-terminator to disable gcc warning
                break;
            case '?':
            error: fprintf(stderr,
                        "Usage:\n"
                        "-d \t number of times a ray can bounce\n"
                        "-s \t number of samples per pixel\n"
                        "-r \t image resolution to be computed\n"
                        "-f \t output file name\n"
                        "-n \t no output(default: 0)\n"
                        "\n"
                        "Example:\n"
                        "%s -d 10 -s 50 -r 720x480 -f tracer.ppm\n",
                        argv[0]);
                exit(EXIT_FAILURE);
                break;
        }
    }

    // Calculating the aspect ratio and creating the camera for the rendering


    if (!no_output)
        fprintf(stderr, "Output file: %s\n", file_name);
    else {
        fprintf(stderr, "No output will be written\n");
    }

    readInput();
    create_random_scene(spheres);

    image_data = static_cast<int*>(malloc(width * height * sizeof(int) * 3));

    for(int i=0; i<THREADS; i++){
	thread_number[i] = i;
	threads[i] = std::thread(kernel, thread_number[i]);
    }
    for(int i=0; i<THREADS; i++){
	threads[i].join();
    }
    
    //Saving the render with PPM format
    if(!no_output) {
        FILE* file;
        if ((file = fopen(file_name, "w")) == NULL )
        {
            perror("fopen");
            exit(EXIT_FAILURE);
        }
        if (fprintf(file, "P3\n%d %d %d\n", width, height, 255) < 0)
        {
            perror("fprintf");
            exit(EXIT_FAILURE);
        }
        for(int y = 0; y < height; y++) {
            for(int x = 0; x < width; x++) {
                int pos = (y * width + x) * 3;
                if (fprintf(file, "%d %d %d\n", image_data[pos] , image_data[pos + 1] , image_data[pos + 2] ) < 0)
                {
                    perror("fprintf");
                    exit(EXIT_FAILURE);
                }
            }
        }
        fclose(file);
    }
    writeOutput(checksum);
    free(image_data);
    return 0;
}
