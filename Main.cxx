#include <iostream>
#include <vector>
#include <set>

#include "SimplRay.hxx"

using namespace RayTracer;

static void WriteToXpm(std::ostream* out, const RasterImage& image);

int main(int argc, char** argv)
{
  std::vector<const ISurface*> surfaces;
  std::vector<Light> lights;

  surfaces.push_back(new Plane
  (
      Vector3F(1000, -100, 100),
      Vector3F(0, 1, -0.2),
      Material ( ColorF(0.5F, 0.5F, 0.5F), 1.0F, 2.0F, 5.0F )
  ));
  surfaces.push_back(new Sphere
  (
      Vector3F(233, 290, 10),
      100,
      Material ( ColorF(1,1,0), 0.5F, 2.0F, 5.0F )
  ));
  surfaces.push_back(new Sphere
  (
      Vector3F(407, 290, -10),
      100,
      Material ( ColorF(0, 1, 1), 0.5F, 2.0F, 5.0F )
  ));
  surfaces.push_back(new Sphere
  (
      Vector3F(320, 140, 0),
      100,
      Material ( ColorF(1, 0, 1), 0.5F, 2.0F, 5.0F )
  ));
  surfaces.push_back(new Sphere
  (
      Vector3F(150, 200, -75),
      10,
      Material ( ColorF(0, 0.5F, 0), 0.5F, 2.0F, 5.0F )
  ));
 

  lights.push_back(Light
  (
       Vector3F(0, 240, -100),
       ColorF(1, 1, 1)
  ));
  lights.push_back(Light
  (
      Vector3F(640, 240, -10000),
      (ColorF(0.6F, 0.7F, 1)).Scale(0.01F)
  ));

  Scene scene;
  scene.ViewWidth = 1280;
  scene.ViewHeight = 1024;
  scene.Mag = 2.0;
  scene.Surfaces = surfaces;
  scene.Lights = lights;
  
  std::cerr << "Generating scene ...." << std::endl;  
  SimplRay simpRay(scene);   
  RasterImage image = simpRay.generate();
  std::cerr << "Writing output ...." << std::endl;
  
  WriteToXpm(&std::cout,image);     
}

static void WriteToXpm(std::ostream* out, const RasterImage& image)
{
  std::set<std::string> colorSet;
  for(int y = 0;y<image.height();++y)
  {
    for(int x = 0;x<image.width();++x)
    {
      colorSet.insert(image.GetPixel(x,y).Scale(255.0).ToRGBHexString());
    }
  }

  (*out) << "/* XPM */n";
  (*out) << "static char * out_xpm[] = {\n";
  (*out) << "\"" << image.width() << " " << image.height() << " " 
         <<  colorSet.size() << " 6\",\n";
  
  for(std::set<std::string>::const_iterator iter = colorSet.begin();
    iter != colorSet.end();
    ++iter)
  {
    (*out) << "\"" << *iter << " c #" << *iter << "\",\n";
  }
  
  for(int y=0;y<image.height();++y)
  {
    if(y!=0) (*out) << ",\n";
    (*out) << "\"";
    for(int x= 0;x<image.width();++x)
    {
      (*out) << image.GetPixel(x,y).Scale(255.0).ToRGBHexString();
    }
    (*out) << "\"";
  }
  
  (*out) << "\n};\n" << std::endl;    
}

