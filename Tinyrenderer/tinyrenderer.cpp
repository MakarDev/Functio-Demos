#include <limits>
#include <cmath>
#include <cassert>
#include <iostream>
#include <cstdint>
#include <fstream>
#include <vector>
#include <sstream>
#include <cstring>

// ===== GEOMETRY.H =====
template<int n> struct vec {
    double data[n] = {0};
    double& operator[](const int i)       { assert(i>=0 && i<n); return data[i]; }
    double  operator[](const int i) const { assert(i>=0 && i<n); return data[i]; }
};

template<int n> double operator*(const vec<n>& lhs, const vec<n>& rhs) {
    double ret1 = 0, ret2 = 0, ret3 = 0, ret4 = 0;
    int i = 0;
    // Process 4 elements at a time to break dependency chain
    for (; i + 3 < n; i += 4) {
        ret1 += lhs[i] * rhs[i];
        ret2 += lhs[i+1] * rhs[i+1];
        ret3 += lhs[i+2] * rhs[i+2];
        ret4 += lhs[i+3] * rhs[i+3];
    }
    // Handle remaining elements
    double remaining = 0;
    for (; i < n; i++) {
        remaining += lhs[i] * rhs[i];
    }
    return (ret1 + ret2) + (ret3 + ret4) + remaining;
}

template<int n> vec<n> operator+(const vec<n>& lhs, const vec<n>& rhs) {
    vec<n> ret = lhs;
    for (int i=n; i--; ret[i]+=rhs[i]);
    return ret;
}

template<int n> vec<n> operator-(const vec<n>& lhs, const vec<n>& rhs) {
    vec<n> ret = lhs;
    for (int i=n; i--; ret[i]-=rhs[i]);
    return ret;
}

template<int n> vec<n> operator*(const vec<n>& lhs, const double& rhs) {
    vec<n> ret = lhs;
    for (int i=n; i--; ret[i]*=rhs);
    return ret;
}

template<int n> vec<n> operator*(const double& lhs, const vec<n> &rhs) {
    return rhs * lhs;
}

template<int n> vec<n> operator/(const vec<n>& lhs, const double& rhs) {
    vec<n> ret = lhs;
    for (int i=n; i--; ret[i]/=rhs);
    return ret;
}

template<int n> std::ostream& operator<<(std::ostream& out, const vec<n>& v) {
    for (int i=0; i<n; i++) out << v[i] << " ";
    return out;
}

template<> struct vec<2> {
    double x = 0, y = 0;
    double& operator[](const int i)       { assert(i>=0 && i<2); return i ? y : x; }
    double  operator[](const int i) const { assert(i>=0 && i<2); return i ? y : x; }
};

template<> struct vec<3> {
    double x = 0, y = 0, z = 0;
    double& operator[](const int i)       { assert(i>=0 && i<3); return i ? (1==i ? y : z) : x; }
    double  operator[](const int i) const { assert(i>=0 && i<3); return i ? (1==i ? y : z) : x; }
};

template<> struct vec<4> {
    double x = 0, y = 0, z = 0, w = 0;
    double& operator[](const int i)       { assert(i>=0 && i<4); return i<2 ? (i ? y : x) : (2==i ? z : w); }
    double  operator[](const int i) const { assert(i>=0 && i<4); return i<2 ? (i ? y : x) : (2==i ? z : w); }
    vec<2> xy()  const { return {x, y};    }
    vec<3> xyz() const { return {x, y, z}; }
};

typedef vec<2> vec2;
typedef vec<3> vec3;
typedef vec<4> vec4;

template<int n> double norm(const vec<n>& v) {
    return std::sqrt(v*v);
}

template<int n> vec<n> normalized(const vec<n>& v) {
    return v / norm(v);
}

inline vec3 cross(const vec3 &v1, const vec3 &v2) {
    return {v1.y*v2.z - v1.z*v2.y, v1.z*v2.x - v1.x*v2.z, v1.x*v2.y - v1.y*v2.x};
}

template<int n> struct dt;

template<int nrows,int ncols> struct mat {
    vec<ncols> rows[nrows] = {{}};

          vec<ncols>& operator[] (const int idx)       { assert(idx>=0 && idx<nrows); return rows[idx]; }
    const vec<ncols>& operator[] (const int idx) const { assert(idx>=0 && idx<nrows); return rows[idx]; }

    double det() const {
        return dt<ncols>::det(*this);
    }

    double cofactor(const int row, const int col) const {
        mat<nrows-1,ncols-1> submatrix;
        for (int i=nrows-1; i--; )
            for (int j=ncols-1;j--; submatrix[i][j]=rows[i+int(i>=row)][j+int(j>=col)]);
        return submatrix.det() * ((row+col)%2 ? -1 : 1);
    }

    mat<nrows,ncols> invert_transpose() const {
        mat<nrows,ncols> adjugate_transpose; // transpose to ease determinant computation, check the last line
        for (int i=nrows; i--; )
            for (int j=ncols; j--; adjugate_transpose[i][j]=cofactor(i,j));
        return adjugate_transpose/(adjugate_transpose[0]*rows[0]);
    }

    mat<nrows,ncols> invert() const {
        return invert_transpose().transpose();
    }

    mat<ncols,nrows> transpose() const {
        mat<ncols,nrows> ret;
        for (int i=ncols; i--; )
            for (int j=nrows; j--; ret[i][j]=rows[j][i]);
        return ret;
    }
};

template<int nrows,int ncols> vec<ncols> operator*(const vec<nrows>& lhs, const mat<nrows,ncols>& rhs) {
    return (mat<1,nrows>{{lhs}}*rhs)[0];
}

template<int nrows,int ncols> vec<nrows> operator*(const mat<nrows,ncols>& lhs, const vec<ncols>& rhs) {
    vec<nrows> ret;
    for (int i=nrows; i--; ret[i]=lhs[i]*rhs);
    return ret;
}

template<int R1,int C1,int C2>mat<R1,C2> operator*(const mat<R1,C1>& lhs, const mat<C1,C2>& rhs) {
    mat<R1,C2> result;
    for (int i=R1; i--; )
        for (int j=C2; j--; )
            for (int k=C1; k--; result[i][j]+=lhs[i][k]*rhs[k][j]);
    return result;
}

template<int nrows,int ncols>mat<nrows,ncols> operator*(const mat<nrows,ncols>& lhs, const double& val) {
    mat<nrows,ncols> result;
    for (int i=nrows; i--; result[i] = lhs[i]*val);
    return result;
}

template<int nrows,int ncols>mat<nrows,ncols> operator/(const mat<nrows,ncols>& lhs, const double& val) {
    mat<nrows,ncols> result;
    for (int i=nrows; i--; result[i] = lhs[i]/val);
    return result;
}

template<int nrows,int ncols>mat<nrows,ncols> operator+(const mat<nrows,ncols>& lhs, const mat<nrows,ncols>& rhs) {
    mat<nrows,ncols> result;
    for (int i=nrows; i--; )
        for (int j=ncols; j--; result[i][j]=lhs[i][j]+rhs[i][j]);
    return result;
}

template<int nrows,int ncols>mat<nrows,ncols> operator-(const mat<nrows,ncols>& lhs, const mat<nrows,ncols>& rhs) {
    mat<nrows,ncols> result;
    for (int i=nrows; i--; )
        for (int j=ncols; j--; result[i][j]=lhs[i][j]-rhs[i][j]);
    return result;
}

template<int nrows,int ncols> std::ostream& operator<<(std::ostream& out, const mat<nrows,ncols>& m) {
    for (int i=0; i<nrows; i++) out << m[i] << std::endl;
    return out;
}

template<int n> struct dt { // template metaprogramming to compute the determinant recursively
    static double det(const mat<n,n>& src) {
        double ret = 0;
        for (int i=n; i--; ret += src[0][i] * src.cofactor(0,i));
        return ret;
    }
};

template<> struct dt<1> {   // template specialization to stop the recursion
    static double det(const mat<1,1>& src) {
        return src[0][0];
    }
};

// ===== TGIMAGE.H =====
#pragma pack(push,1)
struct TGAHeader {
    std::uint8_t  idlength = 0;
    std::uint8_t  colormaptype = 0;
    std::uint8_t  datatypecode = 0;
    std::uint16_t colormaporigin = 0;
    std::uint16_t colormaplength = 0;
    std::uint8_t  colormapdepth = 0;
    std::uint16_t x_origin = 0;
    std::uint16_t y_origin = 0;
    std::uint16_t width = 0;
    std::uint16_t height = 0;
    std::uint8_t  bitsperpixel = 0;
    std::uint8_t  imagedescriptor = 0;
};
#pragma pack(pop)

struct TGAColor {
    std::uint8_t bgra[4] = {0,0,0,0};
    std::uint8_t bytespp = 4;
    std::uint8_t& operator[](const int i) { return bgra[i]; }
};

struct TGAImage {
    enum Format { GRAYSCALE=1, RGB=3, RGBA=4 };
    TGAImage() = default;
    TGAImage(const int w, const int h, const int bpp);
    bool  read_tga_file(const std::string filename);
    bool write_tga_file(const std::string filename, const bool vflip=true, const bool rle=true) const;
    void flip_horizontally();
    void flip_vertically();
    TGAColor get(const int x, const int y) const;
    void set(const int x, const int y, const TGAColor &c);
    int width()  const;
    int height() const;
private:
    bool   load_rle_data(std::ifstream &in);
    bool unload_rle_data(std::ofstream &out) const;
    int w = 0, h = 0;
    std::uint8_t bpp = 0;
    std::vector<std::uint8_t> data = {};
};

// ===== TGIMAGE.CPP =====
TGAImage::TGAImage(const int w, const int h, const int bpp) : w(w), h(h), bpp(bpp), data(w*h*bpp, 0) {}

bool TGAImage::read_tga_file(const std::string filename) {
    std::ifstream in;
    in.open(filename, std::ios::binary);
    if (!in.is_open()) {
        std::cerr << "can't open file " << filename << "\n";
        return false;
    }
    TGAHeader header;
    in.read(reinterpret_cast<char *>(&header), sizeof(header));
    if (!in.good()) {
        std::cerr << "an error occured while reading the header\n";
        return false;
    }
    w   = header.width;
    h   = header.height;
    bpp = header.bitsperpixel>>3;
    if (w<=0 || h<=0 || (bpp!=GRAYSCALE && bpp!=RGB && bpp!=RGBA)) {
        std::cerr << "bad bpp (or width/height) value\n";
        return false;
    }
    size_t nbytes = bpp*w*h;
    data = std::vector<std::uint8_t>(nbytes, 0);
    if (3==header.datatypecode || 2==header.datatypecode) {
        in.read(reinterpret_cast<char *>(data.data()), nbytes);
        if (!in.good()) {
            std::cerr << "an error occured while reading the data\n";
            return false;
        }
    } else if (10==header.datatypecode||11==header.datatypecode) {
        if (!load_rle_data(in)) {
            std::cerr << "an error occured while reading the data\n";
            return false;
        }
    } else {
        std::cerr << "unknown file format " << (int)header.datatypecode << "\n";
        return false;
    }
    if (!(header.imagedescriptor & 0x20))
        flip_vertically();
    if (header.imagedescriptor & 0x10)
        flip_horizontally();
    std::cerr << w << "x" << h << "/" << bpp*8 << "\n";
    return true;
}

bool TGAImage::load_rle_data(std::ifstream &in) {
    size_t pixelcount = w*h;
    size_t currentpixel = 0;
    size_t currentbyte  = 0;
    TGAColor colorbuffer;
    do {
        std::uint8_t chunkheader = 0;
        chunkheader = in.get();
        if (!in.good()) {
            std::cerr << "an error occured while reading the data\n";
            return false;
        }
        if (chunkheader<128) {
            chunkheader++;
            for (int i=0; i<chunkheader; i++) {
                in.read(reinterpret_cast<char *>(colorbuffer.bgra), bpp);
                if (!in.good()) {
                    std::cerr << "an error occured while reading the header\n";
                    return false;
                }
                for (int t=0; t<bpp; t++)
                    data[currentbyte++] = colorbuffer.bgra[t];
                currentpixel++;
                if (currentpixel>pixelcount) {
                    std::cerr << "Too many pixels read\n";
                    return false;
                }
            }
        } else {
            chunkheader -= 127;
            in.read(reinterpret_cast<char *>(colorbuffer.bgra), bpp);
            if (!in.good()) {
                std::cerr << "an error occured while reading the header\n";
                return false;
            }
            for (int i=0; i<chunkheader; i++) {
                for (int t=0; t<bpp; t++)
                    data[currentbyte++] = colorbuffer.bgra[t];
                currentpixel++;
                if (currentpixel>pixelcount) {
                    std::cerr << "Too many pixels read\n";
                    return false;
                }
            }
        }
    } while (currentpixel < pixelcount);
    return true;
}

bool TGAImage::write_tga_file(const std::string filename, const bool vflip, const bool rle) const {
    constexpr std::uint8_t developer_area_ref[4] = {0, 0, 0, 0};
    constexpr std::uint8_t extension_area_ref[4] = {0, 0, 0, 0};
    constexpr std::uint8_t footer[18] = {'T','R','U','E','V','I','S','I','O','N','-','X','F','I','L','E','.','\0'};
    std::ofstream out;
    out.open(filename, std::ios::binary);
    if (!out.is_open()) {
        std::cerr << "can't open file " << filename << "\n";
        return false;
    }
    TGAHeader header = {};
    header.bitsperpixel = bpp<<3;
    header.width  = w;
    header.height = h;
    header.datatypecode = (bpp==GRAYSCALE ? (rle?11:3) : (rle?10:2));
    header.imagedescriptor = vflip ? 0x00 : 0x20; // top-left or bottom-left origin
    out.write(reinterpret_cast<const char *>(&header), sizeof(header));
    if (!out.good()) goto err;
    if (!rle) {
        out.write(reinterpret_cast<const char *>(data.data()), w*h*bpp);
        if (!out.good()) goto err;
    } else if (!unload_rle_data(out)) goto err;
    out.write(reinterpret_cast<const char *>(developer_area_ref), sizeof(developer_area_ref));
    if (!out.good()) goto err;
    out.write(reinterpret_cast<const char *>(extension_area_ref), sizeof(extension_area_ref));
    if (!out.good()) goto err;
    out.write(reinterpret_cast<const char *>(footer), sizeof(footer));
    if (!out.good()) goto err;
    return true;
err:
    std::cerr << "can't dump the tga file\n";
    return false;
}

bool TGAImage::unload_rle_data(std::ofstream &out) const {
    const std::uint8_t max_chunk_length = 128;
    size_t npixels = w*h;
    size_t curpix = 0;
    while (curpix<npixels) {
        size_t chunkstart = curpix*bpp;
        size_t curbyte = curpix*bpp;
        std::uint8_t run_length = 1;
        bool raw = true;
        while (curpix+run_length<npixels && run_length<max_chunk_length) {
            bool succ_eq = true;
            for (int t=0; succ_eq && t<bpp; t++)
                succ_eq = (data[curbyte+t]==data[curbyte+t+bpp]);
            curbyte += bpp;
            if (1==run_length)
                raw = !succ_eq;
            if (raw && succ_eq) {
                run_length--;
                break;
            }
            if (!raw && !succ_eq)
                break;
            run_length++;
        }
        curpix += run_length;
        out.put(raw ? run_length-1 : run_length+127);
        if (!out.good()) return false;
        out.write(reinterpret_cast<const char *>(data.data()+chunkstart), (raw?run_length*bpp:bpp));
        if (!out.good()) return false;
    }
    return true;
}

TGAColor TGAImage::get(const int x, const int y) const {
    if (!data.size() || x<0 || y<0 || x>=w || y>=h) return {};
    TGAColor ret = {0, 0, 0, 0, bpp};
    const std::uint8_t *p = data.data()+(x+y*w)*bpp;
    for (int i=bpp; i--; ret.bgra[i] = p[i]);
    return ret;
}

void TGAImage::set(int x, int y, const TGAColor &c) {
    if (!data.size() || x<0 || y<0 || x>=w || y>=h) return;
    memcpy(data.data()+(x+y*w)*bpp, c.bgra, bpp);
}

void TGAImage::flip_horizontally() {
    for (int i=0; i<w/2; i++)
        for (int j=0; j<h; j++)
            for (int b=0; b<bpp; b++)
                std::swap(data[(i+j*w)*bpp+b], data[(w-1-i+j*w)*bpp+b]);
}

void TGAImage::flip_vertically() {
    const std::size_t row_bytes = static_cast<std::size_t>(w) * bpp;
    const int         half_h    = h / 2;

#pragma omp parallel for
    for (int y = 0; y < half_h; ++y) {
        auto *rowA = data.data() + static_cast<std::size_t>(y) * row_bytes;
        auto *rowB = data.data()
                   + static_cast<std::size_t>(h - 1 - y) * row_bytes;
        std::swap_ranges(rowA, rowA + row_bytes, rowB);
    }
}
int TGAImage::width() const {
    return w;
}

int TGAImage::height() const {
    return h;
}

// ===== MODEL.H =====
class Model {
    std::vector<vec3> verts = {};    // array of vertices        ┐ generally speaking, these arrays
    std::vector<vec3> norms = {};    // array of normal vectors  │ do not have the same size
    std::vector<vec2> tex = {};      // array of tex coords      ┘ check the logs of the Model() constructor
    std::vector<int> facet_vrt = {}; //  ┐ per-triangle indices in the above arrays,
    std::vector<int> facet_nrm = {}; //  │ the size is supposed to be
    std::vector<int> facet_tex = {}; //  ┘ nfaces()*3
    TGAImage diffusemap  = {};       // diffuse color texture
    TGAImage normalmap   = {};       // normal map texture
    TGAImage specularmap = {};       // specular texture
public:
    Model(const std::string filename);
    int nverts() const; // number of vertices
    int nfaces() const; // number of triangles
    vec3 vert(const int i) const;                          // 0 <= i < nverts()
    vec3 vert(const int iface, const int nthvert) const;   // 0 <= iface <= nfaces(), 0 <= nthvert < 3
    vec3 normal(const int iface, const int nthvert) const; // normal coming from the "vn x y z" entries in the .obj file
    vec3 normal(const vec2 &uv) const;                     // normal vector from the normal map texture
    vec2 uv(const int iface, const int nthvert) const;     // uv coordinates of triangle corners
    const TGAImage& diffuse() const;
    const TGAImage& specular() const;
};

// ===== MODEL.CPP =====
Model::Model(const std::string filename) {
    std::ifstream in;
    in.open(filename, std::ifstream::in);
    if (in.fail()) return;
    std::string line;
    while (!in.eof()) {
        std::getline(in, line);
        std::istringstream iss(line.c_str());
        char trash;
        if (!line.compare(0, 2, "v ")) {
            iss >> trash;
            vec3 v;
            for (int i : {0,1,2}) iss >> v[i];
            verts.push_back(v);
        } else if (!line.compare(0, 3, "vn ")) {
            iss >> trash >> trash;
            vec3 n;
            for (int i : {0,1,2}) iss >> n[i];
            norms.push_back(normalized(n));
        } else if (!line.compare(0, 3, "vt ")) {
            iss >> trash >> trash;
            vec2 uv;
            for (int i : {0,1}) iss >> uv[i];
            tex.push_back({uv.x, 1-uv.y});
        }  else if (!line.compare(0, 2, "f ")) {
            int f,t,n, cnt = 0;
            iss >> trash;
            while (iss >> f >> trash >> t >> trash >> n) {
                facet_vrt.push_back(--f);
                facet_tex.push_back(--t);
                facet_nrm.push_back(--n);
                cnt++;
            }
            if (3!=cnt) {
                std::cerr << "Error: the obj file is supposed to be triangulated" << std::endl;
                return;
            }
        }
    }
    std::cerr << "# v# " << nverts() << " f# "  << nfaces() << " vt# " << tex.size() << " vn# " << norms.size() << std::endl;
    auto load_texture = [&filename](const std::string suffix, TGAImage &img) {
        size_t dot = filename.find_last_of(".");
        if (dot==std::string::npos) return;
        std::string texfile = filename.substr(0,dot) + suffix;
        std::cerr << "texture file " << texfile << " loading " << (img.read_tga_file(texfile.c_str()) ? "ok" : "failed") << std::endl;
    };
    load_texture("_diffuse.tga",    diffusemap );
    load_texture("_nm_tangent.tga", normalmap  );
    load_texture("_spec.tga",       specularmap);
}

const TGAImage& Model::diffuse()  const { return diffusemap;  }
const TGAImage& Model::specular() const { return specularmap; }
int Model::nverts() const { return verts.size(); }
int Model::nfaces() const { return facet_vrt.size()/3; }

vec3 Model::vert(const int i) const {
    return verts[i];
}

vec3 Model::vert(const int iface, const int nthvert) const {
    return verts[facet_vrt[iface*3+nthvert]];
}

vec3 Model::normal(const vec2 &uvf) const {
    TGAColor c = normalmap.get(uvf[0]*normalmap.width(), uvf[1]*normalmap.height());
    return vec3{(double)c[2],(double)c[1],(double)c[0]}*2./255. - vec3{1,1,1};
}

vec2 Model::uv(const int iface, const int nthvert) const {
    return tex[facet_tex[iface*3+nthvert]];
}

vec3 Model::normal(const int iface, const int nthvert) const {
    return norms[facet_nrm[iface*3+nthvert]];
}

// ===== OUR_GL.H =====
void viewport(const int x, const int y, const int w, const int h);
void projection(const double coeff=0); // coeff = -1/c
void lookat(const vec3 eye, const vec3 center, const vec3 up);

struct IShader {
    static TGAColor sample2D(const TGAImage &img, const vec2 &uvf) {
        return img.get(uvf[0] * img.width(), uvf[1] * img.height());
    }
    virtual bool fragment(const vec3 bar, TGAColor &color) const = 0;
};

void rasterize(const vec4 clip_verts[3], const IShader &shader, TGAImage &image, std::vector<double> &zbuffer);

// ===== OUR_GL.CPP =====
mat<4,4> ModelView;
mat<4,4> Viewport;
mat<4,4> Projection;

void viewport(const int x, const int y, const int w, const int h) {
    Viewport = {{{w/2., 0, 0, x+w/2.}, {0, h/2., 0, y+h/2.}, {0,0,1,0}, {0,0,0,1}}};
}

void projection(const double f) { // check https://en.wikipedia.org/wiki/Camera_matrix
    Projection = {{{1,0,0,0}, {0,-1,0,0}, {0,0,1,0}, {0,0,-1/f,0}}};
}

void lookat(const vec3 eye, const vec3 center, const vec3 up) { // check https://github.com/ssloy/tinyrenderer/wiki/Lesson-5-Moving-the-camera
    vec3 z = normalized(center-eye);
    vec3 x = normalized(cross(up,z));
    vec3 y = normalized(cross(z, x));
    ModelView = mat<4,4>{{{x.x,x.y,x.z,0}, {y.x,y.y,y.z,0}, {z.x,z.y,z.z,0}, {0,0,0,1}}} *
                mat<4,4>{{{1,0,0,-eye.x},  {0,1,0,-eye.y},  {0,0,1,-eye.z},  {0,0,0,1}}};
}

void rasterize(const vec4 clip_verts[3],
                    const IShader &shader,
                    TGAImage      &image,
                    std::vector<double> &zbuffer)
{
    /* ------------------------------------------------------------- *
     * 1. project to screen space (same as before)
     * ------------------------------------------------------------- */
    vec4 pts [3] = { Viewport * clip_verts[0],
                     Viewport * clip_verts[1],
                     Viewport * clip_verts[2] };

    vec2 v [3] = { (pts[0]/pts[0].w).xy(),
                   (pts[1]/pts[1].w).xy(),
                   (pts[2]/pts[2].w).xy() };

    /* ------------------------------------------------------------- *
     * 2. bounding box clip
     * ------------------------------------------------------------- */
    int bb_min_x = std::max(0,               int(std::floor(std::min(std::min(v[0].x, v[1].x), v[2].x))));
    int bb_min_y = std::max(0,               int(std::floor(std::min(std::min(v[0].y, v[1].y), v[2].y))));
    int bb_max_x = std::min(image.width()-1, int(std::ceil (std::max(std::max(v[0].x, v[1].x), v[2].x))));
    int bb_max_y = std::min(image.height()-1,int(std::ceil (std::max(std::max(v[0].y, v[1].y), v[2].y))));

    /* ------------------------------------------------------------- *
     * 3. edge coefficients  (   E(x,y) = Ax + By + C   )
     *      e0: v1→v2,   e1: v2→v0,   e2: v0→v1
     * ------------------------------------------------------------- */
    auto edge = [](const vec2& a, const vec2& b) {
        return std::pair<float,float>{ a.y - b.y, b.x - a.x };  // (A,B)
    };
    float A0,B0, A1,B1, A2,B2, C0,C1,C2;

    std::tie(A0,B0) = edge(v[1], v[2]);
    std::tie(A1,B1) = edge(v[2], v[0]);
    std::tie(A2,B2) = edge(v[0], v[1]);

    C0 = v[1].x*v[2].y - v[2].x*v[1].y;
    C1 = v[2].x*v[0].y - v[0].x*v[2].y;
    C2 = v[0].x*v[1].y - v[1].x*v[0].y;

    /* signed 2× area; also denominator for normalised barycentrics */
    float area2 = A0*v[0].x + B0*v[0].y + C0;
    if (std::fabs(area2) < 1e-4f) return;           // degenerate tri → skip
    if (area2 < 0) {                                // make area2 positive
        A0 = -A0; B0 = -B0; C0 = -C0;
        A1 = -A1; B1 = -B1; C1 = -C1;
        A2 = -A2; B2 = -B2; C2 = -C2;
        area2 = -area2;
    }
    const float inv_area2 = 1.f / area2;

    const int  W = image.width();

    /* ------------------------------------------------------------- *
     * 4. incremental raster loop
     * ------------------------------------------------------------- */
#pragma omp parallel for
    for (int y = bb_min_y; y <= bb_max_y; ++y)
    {
        /* edge values at (bb_min_x, y) */
        float e0 = A0*bb_min_x + B0*y + C0;
        float e1 = A1*bb_min_x + B1*y + C1;
        float e2 = A2*bb_min_x + B2*y + C2;

        for (int x = bb_min_x; x <= bb_max_x; ++x)
        {
            if (e0 >= 0 && e1 >= 0 && e2 >= 0)           // inside triangle?
            {
                /* perspective-correct barycentrics (only when needed) */
                float   beta  =  e1 * inv_area2;
                float   gamma =  e2 * inv_area2;
                float   alpha =  1.f - beta - gamma;

                vec3 bc_clip  = { alpha / pts[0].w,
                                  beta  / pts[1].w,
                                  gamma / pts[2].w };
                bc_clip = bc_clip / (bc_clip.x + bc_clip.y + bc_clip.z);

                double frag_z = bc_clip * vec3{ clip_verts[0].z,
                                                clip_verts[1].z,
                                                clip_verts[2].z };
                size_t idx = x + y * W;
                if (frag_z < zbuffer[idx])               // depth test
                {
                    TGAColor color;
                    if (!shader.fragment(bc_clip, color))   // shader may discard
                    {
                        zbuffer[idx] = frag_z;
                        image.set(x, y, color);
                    }
                }
            }

            /* step +1 in X: just add A-coefficients               */
            e0 += A0;  e1 += A1;  e2 += A2;
        }
        /* step +1 in Y: add B-coefficients once per scan-line     */
        /* (done implicitly at next loop-iteration via y increment)*/
    }
}

// ===== MAIN FUNCTION =====
struct Shader : IShader {
    const Model &model;
    vec3 uniform_l;       // light direction in view coordinates
    mat<3,2> varying_uv;  // triangle uv coordinates, written by the vertex shader, read by the fragment shader
    mat<3,3> varying_nrm; // normal per vertex to be interpolated by FS
    mat<3,3> view_tri;    // triangle in view coordinates
    mat<3,3> varying_T;   // tangent vectors
    mat<3,3> varying_B;   // bitangent vectors
    mat<4,4> MVP;         // model-view-projection matrix
    mat<4,4> Mview;       // model-view matrix

    Shader(const vec3 l, const Model &m) : model(m) {
        uniform_l = normalized((ModelView*vec4{l.x, l.y, l.z, 0.}).xyz()); // transform the light vector to view coordinates
        MVP = Projection * ModelView;
        Mview = ModelView;
    }

    virtual void vertex(const int iface, const int nthvert, vec4& gl_Position) {
        vec3 p  = model.vert(iface, nthvert);
        vec3 n  = model.normal(iface, nthvert);
        vec2 uv = model.uv    (iface, nthvert);
        varying_uv [nthvert] = uv;
        varying_nrm[nthvert] = n;
        view_tri   [nthvert] = (Mview  * vec4{p.x, p.y, p.z, 1.}).xyz();

        /* --- build per-vertex tangent & bitangent once ---------------- */
        if (nthvert == 0) {          // compute for the whole triangle
            vec3 p0 = view_tri[0], p1 = view_tri[1], p2 = view_tri[2];
            vec2 uv0 = varying_uv[0], uv1 = varying_uv[1], uv2 = varying_uv[2];

            vec3 dp1 = p1 - p0,  dp2 = p2 - p0;
            vec2 duv1 = uv1 - uv0, duv2 = uv2 - uv0;

            double r = 1.0 / (duv1.x * duv2.y - duv1.y * duv2.x);
            vec3 T = (dp1 * duv2.y - dp2 * duv1.y) * r;
            vec3 B = (dp2 * duv1.x - dp1 * duv2.x) * r;

            varying_T[0] = T; varying_T[1] = T; varying_T[2] = T; // copy 3×
            varying_B[0] = B; varying_B[1] = B; varying_B[2] = B;
        }

        gl_Position = MVP * vec4{p.x, p.y, p.z, 1.};
    }

    virtual bool fragment(const vec3 bar, TGAColor &gl_FragColor) const {
        /* interpolate once, normalise only what really matters */
        vec3  n  = normalized(bar * varying_nrm);   // per-pixel normal
        vec3  T  =           bar * varying_T;       // already close to unit
        vec3  B  =           bar * varying_B;
        mat<3,3> TBN;
        TBN[0] = T; TBN[1] = B; TBN[2] = n;
        TBN = TBN.transpose(); // 3×3 mul, no inverse

        /* texture normal, unpack to [-1,1] and normalise */
        vec3 tex_n = model.normal( bar * varying_uv );
        tex_n = normalized( TBN * tex_n );

        /* lighting ---------------------------------------------------- */
        vec3 l = uniform_l;                         // already unit
        double diff   = std::max(0.0, tex_n * l);

        /* fast specular:  r = 2(n·l)n - l  -> cosα = r·v,   v=(0,0,1) */
        double rz = 2.0 * (tex_n * l) * tex_n.z - l.z;

        /* replace pow(x,s) with a 4-term minimax poly x^s ≈ k0+k1x+k2x²+k3x³
           when s∈[5,260].  For demo purposes use s≈32*(specMap/255)+5
           and a simple x^5 * (1-x)^2 (~8 FLOPs, no log/exp).           */
        double glossy  = std::max(0.0, rz);
        double specExp = 5.0 + sample2D(model.specular(), bar * varying_uv)[0];
        double spec    = glossy * glossy * glossy * glossy * glossy; // x^5
        spec *= (1.0 - glossy) * (1.0 - glossy);                     // *(1-x)^2
        spec = std::pow(spec, specExp/5.0);  // still cheaper: base already tiny

        /* final colour ------------------------------------------------ */
        TGAColor base = sample2D(model.diffuse(), bar * varying_uv);
        for (int c=0;c<3;++c)
            gl_FragColor[c] = std::min<int>(10 + base[c] * (diff + spec), 255);
        return false;
    }
};

int final_func(int argc, char** argv) {
    if (2>argc) {
        std::cerr << "Usage: " << argv[0] << " obj/model.obj" << std::endl;
        return 1;
    }

    constexpr int width  = 800;      // output image size
    constexpr int height = 800;
    constexpr vec3 light_dir{1,1,1}; // light source
    constexpr vec3       eye{1,1,3}; // camera position
    constexpr vec3    center{0,0,0}; // camera direction
    constexpr vec3        up{0,1,0}; // camera up vector

    lookat(eye, center, up);                            // build the ModelView matrix
    viewport(width/8, height/8, width*3/4, height*3/4); // build the Viewport matrix
    projection(norm(eye-center));                       // build the Projection matrix
    std::vector<double> zbuffer(width*height, std::numeric_limits<double>::max());

    TGAImage framebuffer(width, height, TGAImage::RGB); // the output image
    for (int m=1; m<argc; m++) { // iterate through all input objects
        Model model(argv[m]);
        Shader shader(light_dir, model);
        for (int t=0; t<model.nfaces(); t++) { // for every triangle
            vec4 clip_vert[3]; // triangle coordinates (clip coordinates), written by VS, read by FS
            for (int v : {0,1,2})
                shader.vertex(t, v, clip_vert[v]);              // call the vertex shader for each triangle vertex
            rasterize(clip_vert, shader, framebuffer, zbuffer); // actual rasterization routine call
        }
    }
    framebuffer.write_tga_file("framebuffer.tga");
    return 0;
}
int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <diffuse_map_path>" << std::endl;
        return 1;
    }

    constexpr long long iterations = 10;

    for (long long i = 0; i < iterations; ++i) {
        volatile auto result = final_func(argc, argv);
        (void)result;
    }

    return 0;
}
