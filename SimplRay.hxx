#ifndef _SimplRay_h_
#define _SimplRay_h_

#include <cmath>
#include <limits>
#include <string>
#include <cstdio>

namespace RayTracer
{

class Constants
{
public:
  static const float E;
  static const float MAXF;
};

const float Constants::E = 2.718281828459045F;
const float Constants::MAXF = std::numeric_limits<float>::max();

class ColorF
{
public:
  ColorF()
  {
      R = 0;
      G = 0;
      B = 0;
  }
  
  ColorF(float r, float g, float b)
  {
      R = r;
      G = g;
      B = b;
  }

  float R;
  float G;
  float B;

  ColorF Clamp(float max) const
  {
      return ColorF(
          std::min(max, R),
          std::min(max, G),
          std::min(max, B));
  }

  ColorF Scale(float x) const
  {
      return ColorF(x*R,x*G,x*B);
  }

  ColorF Filter(const ColorF& f) const
  {
      return ColorF(f.R * R, f.G * G, f.B * B);
  }

  std::string ToRGBHexString() const
  {
    char buf[16];
    sprintf(buf,"%02x%02x%02x",
       std::min(255,(int)R), std::min(255,(int)G), std::min(255,(int)B));
       return buf;
  }

  ColorF Expose(float exposure) const
  {
      return ColorF((float)(1.0 - powf(Constants::E, R * exposure)),
                    (float)(1.0 - powf(Constants::E, G * exposure)),
                    (float)(1.0 - powf(Constants::E, B * exposure)));
  }

  ColorF sRGBEncode() const
  {
      return ColorF(srgbEncode(R), srgbEncode(G), srgbEncode(B));
  }

private:
  static float srgbEncode(float c)
  {
      if (c <= 0.0031308F)
      {
          return 12.92F * c; 
      }
      else
      {
          return (float)(1.055 * powf(c, 0.4166667) - 0.055);
           // Inverse gamma 2.4
      }
  }
};

ColorF operator +(const ColorF& a, const ColorF& b)
{
    return ColorF(a.R + b.R, a.G + b.G, a.B + b.B);
}

class Vector3F
{
public:
  float X;
  float Y;
  float Z;

  Vector3F()
  {
      X = Y = X = 0;
  }

  Vector3F(float x, float y, float z)
  {
      X = x;
      Y = y;
      Z = z;
  }

  operator bool() const { return !(X == 0 && Y == 0 && Z == 0);}
  
  friend Vector3F operator+(const Vector3F& a, const Vector3F& b);
  friend Vector3F operator-(const Vector3F& a, const Vector3F& b);
  friend Vector3F operator*(const Vector3F& a, float r);
  friend float    operator*(const Vector3F& a, const Vector3F& b);
  friend Vector3F operator-(const Vector3F& b);

  Vector3F Dir() const
  {
      if (!(*this)) return Vector3F();

      double mag2 = (*this) * (*this);
      double mag = std::sqrt(mag2);                
      return (*this) * (float)(1.0 / mag);
  }
};

Vector3F operator+(const Vector3F& a, const Vector3F& b)
{
    return Vector3F(a.X + b.X, a.Y + b.Y, a.Z + b.Z);
}
Vector3F operator-(const Vector3F& a, const Vector3F& b)
{
    return Vector3F(a.X - b.X, a.Y - b.Y, a.Z - b.Z);
}
Vector3F operator*(const Vector3F& a, float r)
{
    return Vector3F(r * a.X, r * a.Y, r * a.Z);
}
float operator*(const Vector3F& a, const Vector3F& b)
{
    return a.X * b.X + a.Y * b.Y + a.Z * b.Z;
}
Vector3F operator-(const Vector3F& b)
{
    return Vector3F(- b.X,  - b.Y,  - b.Z);
}

struct Material
{
  Material()
  : Diffuse()
  ,Reflection(0)
  ,SpecValue(0)
  ,SpecPower(0)
  {}
  
  Material(const ColorF& diffuse,float reflection,float specValue,float specPower)
  : Diffuse(diffuse)
  ,Reflection(reflection)
  ,SpecValue(specValue)
  ,SpecPower(specPower)
  {}
  
  ColorF Diffuse;
  float Reflection;
  float SpecValue;
  float SpecPower;
};

struct Light
{
  Light(const Vector3F& position,const ColorF& intensity)
  : Position(position)
  ,Intensity(intensity)
  {}
  
  Vector3F Position;
  ColorF Intensity;
};


class Ray
{
public:
  Ray()
  {
  }
  
  Ray(const Vector3F& origin, const Vector3F& dir)
  : Origin(origin), Dir(dir)
  {}
  
  Vector3F Origin;
  Vector3F Dir;
  Vector3F Point(float t) const { return Origin+Dir*t; }
};


class IntersectResult
{
public:
  IntersectResult()
  : HasResult(false)
  , Ray()
  , Dist(Constants::MAXF)
  , Normal()
  , Material()
  , IP()
  {}
  
  IntersectResult(const Ray& ray,
      float dist,
      const Vector3F& normal,
      const Material& material,
      const Vector3F& iP)
  : HasResult(true)
  , Ray(ray)
  , Dist(dist)
  , Normal(normal)
  , Material(material)
  , IP(iP)
  {}
  
  bool HasResult;
  RayTracer::Ray Ray;
  float Dist;
  Vector3F Normal;
  RayTracer::Material Material;
  Vector3F IP;
};

class ISurface
{
public:
  virtual ~ISurface() {}
  virtual IntersectResult Intersect(const Ray& ray, float max_t) const=0;
  virtual bool CheckIntersect(const Ray& ray, float max_t) const=0;
};
  
class Sphere : public ISurface
{
public:
  Sphere(const Vector3F& centre,float radius,const Material& material)
  : Centre(centre)
  ,Radius(radius)
  ,Material(material)
  {
  }
  
  Vector3F Centre;
  float Radius;
  RayTracer::Material Material;
  
  IntersectResult Intersect(const Ray& ray, float max_t) const
  {
    float new_t = intersect(ray,max_t);
    if(new_t < Constants::MAXF)
    {
      Vector3F ip = ray.Point(new_t);         
      return IntersectResult(ray,new_t,normalAt(ip),Material,ip);
    }
    return IntersectResult();
  }
  
  bool CheckIntersect(const Ray& ray, float max_t) const
  {
    return intersect(ray, max_t)
       < Constants::MAXF;
  }
  
private:  
  float intersect(const Ray& ray, float max_t) const
  {
      Vector3F dist = Centre - ray.Origin;
      float B = ray.Dir * dist;
      float D = (B * B) - (dist * dist) + (Radius * Radius);
      if (D < 0.0F)
      {
          return Constants::MAXF;
      }
      double t0 = B - std::sqrt(D);
      double t1 = B + std::sqrt(D);
      if ((t0 > 0.1f) && (t0 < max_t))
      {
          return (float)t0;
      }
      if ((t1 > 0.1f) && (t1 < max_t))
      {
          return (float)t1;
      }
      return Constants::MAXF;
  }
  
    
  Vector3F normalAt(const Vector3F& p) const
  {
      return (p - Centre).Dir();
  }
};

class Plane : public ISurface
{
public:
  Plane(const Vector3F& point, const Vector3F& normal, const Material& material)
  : Point(point)
  ,Normal(normal.Dir())
  ,Material(material)
  ,Material2(material)
  {
    Material2.Diffuse = ColorF();
  }
  
  
  Vector3F Point;
  Vector3F Normal;
  RayTracer::Material Material;
  RayTracer::Material Material2;
  
  IntersectResult Intersect(const Ray& ray, float max_t) const
  {
    float new_t = intersect(ray,max_t);
    
    if(new_t < Constants::MAXF)
    {
      Vector3F ip = ray.Point(new_t);
      if(((int)(ip.X/50.0) % 2==1)&&((int)(ip.Z/50.0) % 2==1))
      {
        return IntersectResult(ray,new_t,Normal,Material,ip);
      }
      else
      {
       return IntersectResult(ray,new_t,Normal,Material2,ip);
      }
    }
    return IntersectResult();
  }
  
  bool CheckIntersect(const Ray& ray, float max_t) const
  {
    return intersect(ray, max_t)
       < Constants::MAXF;
  }
private:
  float intersect(const Ray& ray, float max_t) const
  {
    float proj = ray.Dir*Normal;
    if(proj>=0.0F)
    {
      return Constants::MAXF;
    }
    
    float new_t = ((Point - ray.Origin) * Normal) / proj;
    if(new_t > 0.1f && new_t < max_t)
    {
      return new_t;
    }
    
    return Constants::MAXF;
  }  
};  

struct Scene
{
  std::vector<const ISurface*> Surfaces;
  std::vector<Light> Lights;
  int ViewWidth;
  int ViewHeight;
  float Mag;
};


class Intersector
{
public:
  static bool CheckIntersects(
    const std::vector<const ISurface*>& surfaces, 
    const Vector3F& origin,
    const Vector3F& to)
  {
      Vector3F dist = (to-origin);
      float max_t = sqrt(dist*dist);
      Ray ray(origin,(to-origin).Dir());
      for(std::vector<const ISurface*>::const_iterator 
           surfIter = surfaces.begin()
          ; surfIter != surfaces.end();
          ++surfIter)
      {       
          const ISurface& s = **surfIter;
          
          if (s.CheckIntersect(ray,max_t))
          {
              return true;
          }
      }          
      return false;
  }

  static IntersectResult GetClosestIntersect(
      const std::vector<const ISurface*>& surfaces, const Ray& ray)
  {
      IntersectResult ir;
      for(std::vector<const ISurface*>::const_iterator 
           surfIter = surfaces.begin()
          ; surfIter != surfaces.end();
          ++surfIter)
      {     
          const ISurface& s = **surfIter;
          IntersectResult new_ir = s.Intersect(ray,ir.Dist);
          if (new_ir.HasResult)
          {
              ir = new_ir;
          }
      }
      return ir;
  }

};

class RasterImage
{
public:
  RasterImage(int w, int h)
  : _w(w)
  , _h(h)
  ,_raster(w*h,ColorF())
  {
  }
  
  void SetPixel(int x, int y, const ColorF& c)
  {
    _raster.at(y*_w+x) = c;
  }
  
  ColorF GetPixel(int x, int y) const 
  {
    return _raster.at(y*_w+x);
  }
    
  int width() const { return _w; }
  int height() const  { return _h; }
private:    
  std::vector<ColorF> _raster;  
  int _w;
  int _h;
};

class SimplRay
{
public:
  SimplRay(Scene scene)
  : m_reflDepth(2)
  , m_superSample(3)
  , m_exposure(-1.0F)
  , m_lambertianCoeff(0.9F)
  , m_blinnPhonegCoeff(1.0F)
  {
      m_scene = scene;      
  }

  RasterImage generate()
  {
      int h = m_scene.ViewHeight;
      int w = m_scene.ViewWidth;

      RasterImage image = RasterImage(w,h);
     
      for (int x = 0; x < w; x++)
      {
          for (int y = 0; y < h; y++)
          {
              image.SetPixel(x, h-y-1,
                  getFinalColor(x,y, m_superSample,1.0F)
                  .Expose(m_exposure)
                  .sRGBEncode()
                  .Clamp(1));
          }
      }

      return image;
  }
private:
  ColorF getFinalColor(float x, float y, int sampleDepth, float w)
  {
      if (sampleDepth <= 0)
      {
          return getFinalColor(x,y);
      }

      int nextDepth = sampleDepth - 1;
      float nw = w / 2;
      return (getFinalColor(x, y, nextDepth, nw) +
             getFinalColor(x + nw, y, nextDepth, nw) +
             getFinalColor(x, y + nw, nextDepth, nw) +
             getFinalColor(x + nw, y + nw, nextDepth, nw)).Scale(0.25F);
  }

  ColorF getFinalColor(float x, float y)
  {
      return getFinalColor(
        Ray(Vector3F(0,0, -2000),
            Vector3F(x / m_scene.Mag,y / m_scene.Mag, 2000).Dir()),
            0);
  }

  ColorF getFinalColor(const Ray& ray, int depth)
  {
      ColorF final(0, 0, 0);

      IntersectResult ir = Intersector::GetClosestIntersect(
        m_scene.Surfaces, ray);

      if (ir.HasResult)            
      {              
          for(std::vector<Light>::const_iterator lightIter = m_scene.Lights.begin()
          ; lightIter != m_scene.Lights.end();
          ++lightIter)
          {     
              const Light& light = *lightIter;
              
              // Check reflected to us
              Ray lightRay(light.Position,(ir.IP - light.Position).Dir());
              
              if (ir.Normal * lightRay.Dir <= 0)
              {                       
                  // Check not in shadow
                  if (!Intersector::CheckIntersects(
                      m_scene.Surfaces,ir.IP, lightRay.Origin))
                  {
                      final = final 
                        //+ getLambertLighting(light, ir, ray, lightRay)
                        + getPhongLighting(light, ir, ray, lightRay);

                  }
              }
          }

          // process reflections
          if (ir.Material.Reflection > 0 && depth < m_reflDepth)
          {
              float reflet = 2.0f * (ray.Dir * ir.Normal);
              ColorF rc = getFinalColor(
                Ray(ir.IP,(ray.Dir - (ir.Normal * reflet)).Dir()), ++depth);

              final = final + rc.Scale(ir.Material.Reflection);
          }           
      }


      return final;
  }

  ColorF getLambertLighting(const Light& light,const IntersectResult& ir,const Ray& viewRay,const Ray& lightRay)
  {
      float lambert = fabsf(lightRay.Dir * ir.Normal) * m_lambertianCoeff;
      return light.Intensity.Scale(lambert).Filter(ir.Material.Diffuse);
  }

  ColorF getBlinnPhongLighting(const Light& light,const IntersectResult& ir,const Ray& viewRay,const Ray& lightRay)
  {
      Vector3F blinnDir = ((-lightRay.Dir) - viewRay.Dir).Dir();
      float blinnTerm = std::max(blinnDir * ir.Normal, 0.0F);
      blinnTerm = ir.Material.SpecValue * 
          (float)powf(blinnTerm, ir.Material.SpecPower) * 
          m_blinnPhonegCoeff;

      return light.Intensity.Scale(blinnTerm).Filter(ir.Material.Diffuse);
  }

  ColorF getPhongLighting(const Light& light,const IntersectResult& ir,const Ray& viewRay,const Ray& lightRay)
  {
      float reflet = 2.0F * ((-lightRay.Dir) * ir.Normal);
      Vector3F phongDir = (-lightRay.Dir) - (ir.Normal * reflet);
      float phongTerm = std::max(phongDir * viewRay.Dir, 0.0f);
      phongTerm = ir.Material.SpecValue * (float)powf(phongTerm, ir.Material.SpecPower) 
          * m_blinnPhonegCoeff;

      return getLambertLighting(light, ir, viewRay, lightRay).Scale(1.0F + phongTerm);
  }

  Scene m_scene;
  int m_reflDepth;
  int m_superSample;
  float m_exposure;
  float m_lambertianCoeff;
  float m_blinnPhonegCoeff;
};

}

#endif
