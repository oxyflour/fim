#include <clipper.hpp>

#ifndef SVG_H
#define SVG_H

#include <string>

namespace clipper {

using namespace std;
using namespace ClipperLib;

/*
class SVGBuilder
{
  static string ColorToHtml(unsigned clr)
  {
    stringstream ss;
    ss << '#' << hex << std::setfill('0') << setw(6) << (clr & 0xFFFFFF);
    return ss.str();
  }
  //------------------------------------------------------------------------------

  static float GetAlphaAsFrac(unsigned clr)
  {
    return ((float)(clr >> 24) / 255);
  }
  //------------------------------------------------------------------------------

  class StyleInfo
  {
  public:
  PolyFillType pft;
  unsigned brushClr;
  unsigned penClr;
  double penWidth;
  bool showCoords;

  StyleInfo()
  {
    pft = pftNonZero;
    brushClr = 0xFFFFFFCC;
    penClr = 0xFF000000;
    penWidth = 0.8;
    showCoords = false;
  }
  };

  class PolyInfo
  {
    public:
      Paths paths;
    StyleInfo si;

      PolyInfo(Paths paths, StyleInfo style)
      {
          this->paths = paths;
          this->si = style;
      }
  };

  typedef std::vector<PolyInfo> PolyInfoList;

private:
  PolyInfoList polyInfos;
  static const std::string svg_xml_start[];
  static const std::string poly_end[];

public:
  StyleInfo style;

  void AddPaths(Paths& poly)
  {
    if (poly.size() == 0) return;
    polyInfos.push_back(PolyInfo(poly, style));
  }

  bool SaveToFile(const string& filename, double scale = 1.0, int margin = 10)
  {
    //calculate the bounding rect ...
    PolyInfoList::size_type i = 0;
    Paths::size_type j;
    while (i < polyInfos.size())
    {
      j = 0;
      while (j < polyInfos[i].paths.size() &&
        polyInfos[i].paths[j].size() == 0) j++;
      if (j < polyInfos[i].paths.size()) break;
      i++;
    }
    if (i == polyInfos.size()) return false;

    IntRect rec;
    rec.left = polyInfos[i].paths[j][0].X;
    rec.right = rec.left;
    rec.top = polyInfos[i].paths[j][0].Y;
    rec.bottom = rec.top;
    for ( ; i < polyInfos.size(); ++i)
      for (Paths::size_type j = 0; j < polyInfos[i].paths.size(); ++j)
        for (Path::size_type k = 0; k < polyInfos[i].paths[j].size(); ++k)
        {
          IntPoint ip = polyInfos[i].paths[j][k];
          if (ip.X < rec.left) rec.left = ip.X;
          else if (ip.X > rec.right) rec.right = ip.X;
          if (ip.Y < rec.top) rec.top = ip.Y;
          else if (ip.Y > rec.bottom) rec.bottom = ip.Y;
        }

    if (scale == 0) scale = 1.0;
    if (margin < 0) margin = 0;
    rec.left = (cInt)((double)rec.left * scale);
    rec.top = (cInt)((double)rec.top * scale);
    rec.right = (cInt)((double)rec.right * scale);
    rec.bottom = (cInt)((double)rec.bottom * scale);
    cInt offsetX = -rec.left + margin;
    cInt offsetY = -rec.top + margin;

    ofstream file;
    file.open(filename);
    if (!file.is_open()) return false;
    file.setf(ios::fixed);
    file.precision(0);
    file << svg_xml_start[0] <<
      ((rec.right - rec.left) + margin*2) << "px" << svg_xml_start[1] <<
      ((rec.bottom - rec.top) + margin*2) << "px" << svg_xml_start[2] <<
      ((rec.right - rec.left) + margin*2) << " " <<
      ((rec.bottom - rec.top) + margin*2) << svg_xml_start[3];
    setlocale(LC_NUMERIC, "C");
    file.precision(2);

    for (PolyInfoList::size_type i = 0; i < polyInfos.size(); ++i)
  {
      file << " <path d=\"";
    for (Paths::size_type j = 0; j < polyInfos[i].paths.size(); ++j)
      {
        if (polyInfos[i].paths[j].size() < 3) continue;
        file << " M " << ((double)polyInfos[i].paths[j][0].X * scale + offsetX) <<
          " " << ((double)polyInfos[i].paths[j][0].Y * scale + offsetY);
        for (Path::size_type k = 1; k < polyInfos[i].paths[j].size(); ++k)
        {
          IntPoint ip = polyInfos[i].paths[j][k];
          double x = (double)ip.X * scale;
          double y = (double)ip.Y * scale;
          file << " L " << (x + offsetX) << " " << (y + offsetY);
        }
        file << " z";
    }
      file << poly_end[0] << ColorToHtml(polyInfos[i].si.brushClr) <<
    poly_end[1] << GetAlphaAsFrac(polyInfos[i].si.brushClr) <<
        poly_end[2] <<
        (polyInfos[i].si.pft == pftEvenOdd ? "evenodd" : "nonzero") <<
        poly_end[3] << ColorToHtml(polyInfos[i].si.penClr) <<
    poly_end[4] << GetAlphaAsFrac(polyInfos[i].si.penClr) <<
        poly_end[5] << polyInfos[i].si.penWidth << poly_end[6];

        if (polyInfos[i].si.showCoords)
        {
      file << "<g font-family=\"Verdana\" font-size=\"11\" fill=\"black\">\n\n";
      for (Paths::size_type j = 0; j < polyInfos[i].paths.size(); ++j)
      {
        if (polyInfos[i].paths[j].size() < 3) continue;
        for (Path::size_type k = 0; k < polyInfos[i].paths[j].size(); ++k)
        {
          IntPoint ip = polyInfos[i].paths[j][k];
          file << "<text x=\"" << (int)(ip.X * scale + offsetX) <<
          "\" y=\"" << (int)(ip.Y * scale + offsetY) << "\">" <<
          ip.X << "," << ip.Y << "</text>\n";
          file << "\n";
        }
      }
      file << "</g>\n";
        }
    }
    file << "</svg>\n";
    file.close();
    setlocale(LC_NUMERIC, "");
    return true;
  }
}; //SVGBuilder
*/

}

#endif
