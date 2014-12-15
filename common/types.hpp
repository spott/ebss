#pragma once

#include<petsc.h>

//stl:
#include<iostream>
#include<vector>
#include<array>

#include<boost/serialization/access.hpp>
#include<boost/serialization/complex.hpp>

template <typename iterator>
struct Range{
    iterator begin;
    iterator end;
};

struct kBasisID{
    PetscReal k;
    PetscInt l;
};

std::ostream& operator<<( std::ostream& out, const kBasisID a) {
    out << a.k << ", " << a.l;
    return out;
}

bool operator==(const kBasisID &a, const kBasisID &b)
{
    if (a.l != b.l || a.k != b.k )
        return false;
    else
        return true;
}

struct BasisID {
    // j = 1 -> 1/2;
    // j = 3 -> 3/2;
    // j != 1,3 -> 0;
    PetscInt n,l,j,m;
    PetscScalar e;
    bool operator<(const BasisID b) const
    {
        if (this->l < b.l)
            return true;
        else if (this->l == b.l && this->n < b.n)
            return true;
        else if (this->l == b.l && this->n < b.n && this->j == b.j )
            return true;
        else if (this->l == b.l && this->n == b.n && this->j == b.j && this->m < b.m)
            return true;
        else
            return false;
    }

//private:
    //friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int /*version*/)
    {
        ar & n;
        ar & l;
        ar & j;
        ar & m;
        ar & e;
    }
};

bool operator==(const BasisID &a, const BasisID &b)
{
    if (a.l != b.l || a.n != b.n || a.j != b.j || a.e != b.e || a.m != b.m )
        return false;
    else
        return true;
}
bool operator!=(const BasisID &a, const BasisID &b)
{
    if (a.l != b.l || a.n != b.n || a.j != b.j || a.e != b.e || a.m != b.m )
        return true;
    else
        return false;
}
std::istream& operator>>(std::istream &in, BasisID &b)     //input
{
    PetscReal er, ei;
    in >> b.n >> b.l >> b.j >> b.m >> er >> ei;
    b.e = std::complex<double>(er,ei);
    return in;
}
std::ostream& operator<<(std::ostream &out, const BasisID &b)     //output
{
    out << b.n << ", " << b.l << ", " << b.j << ", " << b.m << ", " << b.e.real();
    return out;
}

template<typename T, size_t N>
std::ostream& operator<<(std::ostream &out, const std::array<T, N> &b) //output an array
{
    for (size_t i = 0; i < b.size(); i++)
    {
        if (i != 0)
            out << ", " << b[i];
        else
            out << b[i];
    }
    return out;
}

template< typename T, size_t N>
bool operator<( const std::array<T, N> &a, const std::array<T, N> &b )
{
    for (size_t i = 0; i < N; ++i)
    {
        if (a[i] < b[i])
            return true;
        else if (a[i] > b[i])
            return false;
    }
    return false;
}

namespace std {
template<>
struct hash< BasisID > {
public:
    size_t operator()(const BasisID &a) const 
    {
        return a.n + 10000 * a.l + 10000000 * a.j + 100000000 * a.m;
    }
};

template<>
struct hash< kBasisID > {
public:
    size_t operator()(const kBasisID &a) const 
    {
        return 10000000 * a.l + 100000 * a.k;
    }
};
}

template<class T, class A=std::allocator<T> >
struct Matrix {
  typedef T value_type;
  typedef std::vector<value_type, A> Container;

  Matrix() : _b(0) {}
  Matrix(int a, int b, value_type const& initial=value_type())
  : _b(0)
  {
    resize(a, b, initial);
  }
  Matrix(Matrix const& other)
  : _data(other._data), _b(other._b)
  {}

  Matrix& operator=(Matrix copy) {
    swap(*this, copy);
    return *this;
  }

  bool empty() const { return _data.empty(); }
  void clear() { _data.clear(); _b = 0; }

  int dim_x() const { return _b ? _data.size() / _b : 0; }
  int dim_y() const { return _b; }

  value_type* operator[](int a) {
    return &_data[a * _b];
  }
  value_type const* operator[](int a) const {
    return &_data[a * _b];
  }

  value_type& operator()(int x, int y) {
      return _data[x * _b + y];
  }

  void resize(int a, int b, value_type const& initial=value_type()) {
    if (a == 0) {
      b = 0;
    }
    _data.resize(a * b, initial);
    _b = b;
  }

  const Container& data() const { return _data; }

  friend void swap(Matrix& a, Matrix& b) {
    using std::swap;
    swap(a._data, b._data);
    swap(a._b,    b._b);
  }

  template<class Stream>
  friend Stream& operator<<(Stream& s, Matrix const& value) {
    s << "<Matrix at " << &value << " dimensions "
      << value.dim_a() << 'x' << value.dim_b();
    if (!value.empty()) {
      bool first = true;
      typedef typename Container::const_iterator Iter;
      Iter i = value._data.begin(), end = value._data.end();
      while (i != end) {
        s << (first ? " [[" : "], [");
        first = false;
        s << *i;
        ++i;
        for (int b = value._b - 1; b; --b) {
          s << ", " << *i;
          ++i;
        }
      }
      s << "]]";
    }
    s << '>';
    return s;
  }

private:
  Container _data;
  int _b;
};

namespace units {
    //for converting between "real" and "atomic" units:

    class Meter {
        public:
        constexpr Meter( double m ) : v(m) {};

        double value() { return v; };

        private: double v;
    };


    class eV {
        public: 
        constexpr eV( double m ) : v(m) {};

        double value() { return v; };

        private: double v;
    };


    template<typename T>
    double toEnergy( T ) {};

    template<>
    double toEnergy< Meter > ( Meter from ) { return 4.556335e-8 / from.value(); };

    template<>
    double toEnergy< eV > ( eV from ) { return from.value() * 0.03674932; };


}
inline constexpr units::Meter operator"" _nm( long double d ) { return units::Meter( d * 1e-9 ); };
inline constexpr units::Meter operator"" _cm( long double d ) { return units::Meter( d * 1e-2 ); };

inline constexpr units::eV operator"" _eV( long double d ) { return units::eV( d ); };

