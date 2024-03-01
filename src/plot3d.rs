use crate::DrawResult;
use geo::{
    Densify, EuclideanDistance, InteriorPoint, Intersects, LineInterpolatePoint, LineIntersection,
};
use plotters::prelude::*;
use plotters_canvas::CanvasBackend;
use std::collections::BTreeMap;
use std::ops::{Add, Neg};
use std::sync::OnceLock;
use web_sys::HtmlCanvasElement;

static SURFACE_MAP: OnceLock<BTreeMap<(u64, u64), (f64, f64)>> = OnceLock::new();
pub fn draw(canvas: HtmlCanvasElement, pitch: f64, yaw: f64, mvx: f64, mvz: f64) -> DrawResult<()> {
    let area = CanvasBackend::with_canvas_object(canvas)
        .unwrap()
        .into_drawing_area();
    area.fill(&WHITE)?;

    let x_axis = (-3.0..3.0).step(0.1);
    let z_axis = (-3.0..3.0).step(0.1);

    let mut chart =
        ChartBuilder::on(&area).build_cartesian_3d(x_axis.clone(), -3.0..3.0, z_axis.clone())?;

    chart.with_projection(|mut pb| {
        pb.yaw = yaw;
        pb.pitch = pitch;
        pb.scale = 0.7;
        pb.into_matrix()
    });

    chart.configure_axes().draw()?;

    // chart.draw_series(
    //     SurfaceSeries::xoz(x_axis.values(), z_axis.values(), |x: f64, z: f64| {
    //         (x * x + z * z).cos()
    //     })
    //         .style(&BLUE.mix(0.2)),
    // )?;

    // chart.draw_series(LineSeries::new(
    //     (-100..100)
    //         .map(|y| y as f64 / 40.0)
    //         .map(|y| ((y * 10.0).sin(), y, (y * 10.0).cos())),
    //     &BLACK,
    // ))?;

    let point_u = (-2. + mvx, 0. + mvz);
    let point_v = (-1. + mvx, 3.0_f64.sqrt().neg() + mvz);
    let point_w = (1. + mvx, 3.0_f64.sqrt().neg() + mvz);
    let point_x = (2. + mvx, 0. + mvz);
    let point_y = (1. + mvx, 3.0_f64.sqrt() + mvz);
    let point_z = (-1. + mvx, 3.0_f64.sqrt() + mvz);
    let points = [point_u, point_v, point_w, point_x, point_y, point_z];
    let points2 = [point_v, point_w, point_x, point_y, point_z, point_u];
    let lines: Vec<_> = points
        .iter()
        .zip(points2.iter())
        .map(|(&a, &b)| geo::Line::new(a, b))
        .collect();
    let mut line_points = points
        .iter()
        .zip(points2.iter())
        .map(|(&a, &b)| new_line(a, b));
    for line in line_points.clone() {
        chart.draw_series(LineSeries::new(line, &BLACK))?;
    }

    let x_axis = (-1.0..1.0).step(0.1);
    let z_axis = (-1.0..1.0).step(0.1);

    // let mut low_surface: Vec<(f64, f64, f64)> = vec![];
    // let mut high_surface: Vec<(f64, f64, f64)> = vec![];
    // for x in x_axis.values() {
    //     for z in z_axis.values() {
    //         let (min, max) = cacl_min_max_height(&lines, x, z);
    //         low_surface.push((x, min, z));
    //         high_surface.push((x, max, z));
    //     }
    // }
    //
    // chart.draw_series(LineSeries::new(low_surface, &BLUE))?;
    //
    // chart.draw_series(LineSeries::new(high_surface, &BLUE))?;

    let mut surface_map = BTreeMap::new();
    for x in x_axis.values() {
        for z in z_axis.values() {
            let (min, max) = cacl_min_max_height(&lines, x, z);
            let (x, z): (u64, u64) = unsafe { std::mem::transmute((x, z)) };
            surface_map.insert((x, z), (min, max));
        }
    }

    let low_surface = SurfaceSeries::xoz(x_axis.values(), z_axis.values(), |x: f64, z: f64| {
        let (x_u64, z_u64): (u64, u64) = unsafe { std::mem::transmute((x, z)) };
        let y = surface_map.get(&(x_u64, z_u64)).unwrap().0;
        y
    });
    let high_surface = SurfaceSeries::xoz(x_axis.values(), z_axis.values(), |x: f64, z: f64| {
        let (x_u64, z_u64): (u64, u64) = unsafe { std::mem::transmute((x, z)) };
        let y = surface_map.get(&(x_u64, z_u64)).unwrap().1;
        y
    });
    chart.draw_series(low_surface.style(&BLUE.mix(0.2)))?;
    chart.draw_series(high_surface.style(&BLUE.mix(0.2)))?;
    Ok(())
}

fn cacl_min_max_height(lines: &[geo::Line], x: f64, z: f64) -> (f64, f64) {
    let (mut min, mut max) = (4_f64, 0_f64);
    for i in 0..100 {
        let angle = i as f64 / std::f64::consts::PI;
        let h = cacl_height(angle, lines, x, z);
        min = h.min(min);
        max = h.max(max);
    }
    (min, max)
}

fn cacl_height(angle: f64, lines: &[geo::Line], x: f64, z: f64) -> f64 {
    const LEN: f64 = 4.;
    let point_a = geo::Coord::from((LEN * angle.cos(), LEN * angle.sin()));
    let point_a1 = -point_a;
    let point_p = geo::Coord::from((x, z));
    let point_b = point_a + point_p;
    let point_b1 = point_a1 + point_p;
    let line_pb = geo::Line::new(point_p, point_b);
    let mut point_s = None;
    for &line in lines {
        if let Some(LineIntersection::SinglePoint {
            intersection,
            is_proper,
        }) = geo::line_intersection::line_intersection(line, line_pb)
        {
            point_s = Some(intersection);
            break;
        }
    }
    let point_s = point_s.unwrap();
    let point_c = geo::Coord::from((point_a.y, -point_a.x));
    let point_c1 = geo::Coord::from((-point_a.y, point_a.x));
    let line_bb1 = geo::Line::new(point_b, point_b1);
    let line_cc1 = geo::Line::new(point_c, point_c1);
    let point_q = geo::line_intersection::line_intersection(line_bb1, line_cc1).unwrap();
    let LineIntersection::SinglePoint {
        intersection: point_q,
        is_proper,
    } = point_q
    else {
        panic!("")
    };
    let len_qs = point_q.euclidean_distance(&point_s);
    4. - len_qs
}

fn new_line(a: (f64, f64), b: (f64, f64)) -> Vec<(f64, f64, f64)> {
    let line = geo::Line::new(a, b);
    let line_points = (0..=100_u64)
        .map(|v| {
            let point = line.line_interpolate_point(v as f64 / 100.).unwrap();
            (point.x(), 0_f64, point.y())
        })
        .collect();
    line_points
}

#[cfg(test)]
mod test {
    #[test]
    fn test() {}
}
