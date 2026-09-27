//! 2D points and vectors.

use core::ops::{Add, AddAssign, Div, Mul, Neg, Sub, SubAssign};

use crate::sqrt;

/// A 2D point, also used as a vector.
///
/// Deliberately a plain struct: layout and field access are part of the public API, and `f64`
/// semantics are the crate's arithmetic (see `docs/SPEC.md` section 3).
#[derive(Clone, Copy, Debug, Default, PartialEq)]
#[repr(C)]
pub struct Point2D {
    /// `x` coordinate.
    pub x: f64,
    /// `y` coordinate.
    pub y: f64,
}

impl Point2D {
    /// The origin.
    pub const ZERO: Self = Self { x: 0.0, y: 0.0 };

    /// Constructs a point.
    #[inline]
    pub const fn new(x: f64, y: f64) -> Self {
        Self { x, y }
    }

    /// Dot product.
    #[inline]
    pub fn dot(self, other: Self) -> f64 {
        self.x * other.x + self.y * other.y
    }

    /// The scalar cross product (z component of the 3D cross product).
    #[inline]
    pub fn cross(self, other: Self) -> f64 {
        self.x * other.y - self.y * other.x
    }

    /// Squared Euclidean length. Cheaper than [`Self::length`] and enough for comparisons.
    #[inline]
    pub fn norm_squared(self) -> f64 {
        self.dot(self)
    }

    /// Euclidean length.
    #[inline]
    pub fn length(self) -> f64 {
        sqrt(self.norm_squared())
    }

    /// Euclidean distance between two points.
    #[inline]
    pub fn distance(self, other: Self) -> f64 {
        (self - other).length()
    }

    /// Squared Euclidean distance between two points.
    #[inline]
    pub fn distance_squared(self, other: Self) -> f64 {
        (self - other).norm_squared()
    }

    /// Returns the unit vector in this direction, or `None` if this is the zero vector.
    #[inline]
    pub fn normalize(self) -> Option<Self> {
        let len = self.length();
        if len == 0.0 { None } else { Some(self / len) }
    }

    /// Whether both components are finite.
    #[inline]
    pub fn is_finite(self) -> bool {
        self.x.is_finite() && self.y.is_finite()
    }

    /// Whether this is exactly the zero vector (bitwise).
    #[inline]
    pub fn is_zero(self) -> bool {
        self.x == 0.0 && self.y == 0.0
    }

    /// Left-hand perpendicular, `(-y, x)`.
    #[inline]
    pub fn perp(self) -> Self {
        Self::new(-self.y, self.x)
    }

    /// The counter-clockwise perpendicular of `self - other`, unit length when `self != other`.
    #[inline]
    pub fn left_normal_from(self, other: Self) -> Self {
        (self - other).perp().normalize().unwrap_or(Self::ZERO)
    }

    /// Linear interpolation; `t` is not clamped.
    #[inline]
    pub fn lerp(self, other: Self, t: f64) -> Self {
        self + (other - self) * t
    }
}

impl Add for Point2D {
    type Output = Self;
    #[inline]
    fn add(self, rhs: Self) -> Self {
        Self::new(self.x + rhs.x, self.y + rhs.y)
    }
}

impl AddAssign for Point2D {
    #[inline]
    fn add_assign(&mut self, rhs: Self) {
        *self = *self + rhs;
    }
}

impl Sub for Point2D {
    type Output = Self;
    #[inline]
    fn sub(self, rhs: Self) -> Self {
        Self::new(self.x - rhs.x, self.y - rhs.y)
    }
}

impl SubAssign for Point2D {
    #[inline]
    fn sub_assign(&mut self, rhs: Self) {
        *self = *self - rhs;
    }
}

impl Neg for Point2D {
    type Output = Self;
    #[inline]
    fn neg(self) -> Self {
        Self::new(-self.x, -self.y)
    }
}

impl Mul<f64> for Point2D {
    type Output = Self;
    #[inline]
    fn mul(self, rhs: f64) -> Self {
        Self::new(self.x * rhs, self.y * rhs)
    }
}

impl Mul<Point2D> for f64 {
    type Output = Point2D;
    #[inline]
    fn mul(self, rhs: Point2D) -> Point2D {
        rhs * self
    }
}

impl Div<f64> for Point2D {
    type Output = Self;
    #[inline]
    fn div(self, rhs: f64) -> Self {
        Self::new(self.x / rhs, self.y / rhs)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn dot_and_cross() {
        let a = Point2D::new(1.0, 2.0);
        let b = Point2D::new(3.0, 4.0);
        assert_eq!(a.dot(b), 11.0);
        assert_eq!(a.cross(b), -2.0);
    }

    #[test]
    fn length_and_normalize() {
        let a = Point2D::new(3.0, 4.0);
        assert_eq!(a.length(), 5.0);
        assert_eq!(a.normalize(), Some(Point2D::new(0.6, 0.8)));
        assert_eq!(Point2D::ZERO.normalize(), None);
    }

    #[test]
    fn perp_is_left_normal() {
        let a = Point2D::new(1.0, 0.0);
        assert_eq!(a.perp(), Point2D::new(0.0, 1.0));
        assert_eq!(a.left_normal_from(Point2D::ZERO), Point2D::new(0.0, 1.0));
        assert_eq!(a.left_normal_from(a), Point2D::ZERO);
    }

    #[test]
    fn arithmetic_traits() {
        let a = Point2D::new(1.0, 2.0);
        let mut b = Point2D::new(3.0, 4.0);
        assert_eq!(a + b, Point2D::new(4.0, 6.0));
        assert_eq!(b - a, Point2D::new(2.0, 2.0));
        assert_eq!(-a, Point2D::new(-1.0, -2.0));
        assert_eq!(2.0 * a, Point2D::new(2.0, 4.0));
        assert_eq!(a * 2.0, Point2D::new(2.0, 4.0));
        assert_eq!(b / 2.0, Point2D::new(1.5, 2.0));
        b += a;
        assert_eq!(b, Point2D::new(4.0, 6.0));
        b -= a;
        assert_eq!(b, Point2D::new(3.0, 4.0));
    }

    #[test]
    fn finite_and_zero() {
        assert!(Point2D::new(1.0, 2.0).is_finite());
        assert!(!Point2D::new(f64::NAN, 2.0).is_finite());
        assert!(!Point2D::new(f64::INFINITY, 2.0).is_finite());
        assert!(Point2D::ZERO.is_zero());
    }
}
