#ifndef _POINT_H_
#define _POINT_H_ 1

#include <complex>

class Point
{
private:
	double m_x;
	double m_y;

public:
	Point() : m_x(0.0), m_y(0.0) {}
	Point(const double x, const double y) : m_x(x), m_y(y) {}
	~Point() {}

	// accessors
	double& x() { return m_x; }
	double& y() { return m_y; }

	const double& x() const { return m_x; }
	const double& y() const { return m_y; }

  // operator
  bool operator==(const Point& p)const
  {
    return p.x() == m_x && p.y() == m_y;
  };

	// compute inverse distance to q, with power p
	double inverse_distance(const Point& q, const double power) const
	{
		const double d = std::sqrt((q.x() - m_x) * (q.x() - m_x) + 
			                       (q.y() - m_y) * (q.y() - m_y));
		if (d != 0.0)
			return 1.0 / std::pow(d, power);
		else
			return 1.0;
	}

	// return true when p coinces with "this", false otherwise
	bool coincide(const Point& p) const
	{
		return p.x() == m_x && p.y() == m_y;
	}
};

#endif // _POINT_H_

