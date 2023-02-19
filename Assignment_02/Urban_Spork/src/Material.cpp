#include "Material.h"
#include "Primitive.h"

bool Lambert::shade(const Ray& ray, const Intersection& data, Color& attenuation, Ray& scatter) {
    const vec3 target = data.p + data.normal + randomUnitSphere();
    scatter = Ray{data.p, (target - data.p).normalized()};
    attenuation = albedo;
    return true;
}

bool Metal::shade(const Ray& ray, const Intersection& data, Color& attenuation, Ray& scatter) {
    const vec3 reflected =
        (reflect(ray.dir, data.normal).normalized() + fuzz * randomUnitSphere()).normalized();
    scatter = Ray{data.p, reflected};
    attenuation = albedo;
    return dot(scatter.dir, data.normal) > 0.f;
}
