use crate::Primitive;
use crate::core::renderer::Quad;
use crate::core::{
    Background, Color, Gradient, Rectangle, Size, Transformation, Vector,
};
use crate::graphics::{Image, Text};
use crate::text;

use vello_cpu::kurbo::Shape as _;

#[derive(Debug)]
pub struct Engine {
    text_pipeline: text::Pipeline,

    #[cfg(feature = "image")]
    pub(crate) raster_pipeline: crate::raster::Pipeline,
    #[cfg(feature = "svg")]
    pub(crate) vector_pipeline: crate::vector::Pipeline,
}

impl Engine {
    pub fn new() -> Self {
        Self {
            text_pipeline: text::Pipeline::new(),
            #[cfg(feature = "image")]
            raster_pipeline: crate::raster::Pipeline::new(),
            #[cfg(feature = "svg")]
            vector_pipeline: crate::vector::Pipeline::new(),
        }
    }

    pub fn draw_quad(
        &mut self,
        quad: &Quad,
        background: &Background,
        transformation: Transformation,
        pixels: &mut vello_cpu::Pixmap,
        render_context: &mut vello_cpu::RenderContext,
        clip_bounds: Rectangle,
    ) {
        debug_assert!(
            quad.bounds.width.is_normal(),
            "Quad with non-normal width!"
        );
        debug_assert!(
            quad.bounds.height.is_normal(),
            "Quad with non-normal height!"
        );

        let physical_bounds = quad.bounds * transformation;

        if !clip_bounds.intersects(&physical_bounds) {
            return;
        }

        // let clip_mask = (!physical_bounds.is_within(&clip_bounds))
        //     .then_some(clip_mask as &_);

        let transform = into_transform(transformation);

        // Make sure the border radius is not larger than the bounds
        let border_width = quad
            .border
            .width
            .min(quad.bounds.width / 2.0)
            .min(quad.bounds.height / 2.0);

        let mut fill_border_radius = <[f32; 4]>::from(quad.border.radius);

        for radius in &mut fill_border_radius {
            *radius = (*radius)
                .min(quad.bounds.width / 2.0)
                .min(quad.bounds.height / 2.0);
        }

        let path = rounded_rectangle(quad.bounds, fill_border_radius);

        let paint = match background {
            Background::Color(color) => {
                vello_cpu::PaintType::Solid(into_color(*color))
            }
            Background::Gradient(Gradient::Linear(linear)) => {
                let (start, end) = linear.angle.to_distance(&quad.bounds);

                let stops: Vec<vello_cpu::peniko::ColorStop> = linear
                    .stops
                    .into_iter()
                    .flatten()
                    .map(|stop| (stop.offset, into_color(stop.color)).into())
                    .collect();

                vello_cpu::PaintType::Gradient(
                    vello_cpu::peniko::Gradient::new_linear(
                        (start.x, start.y),
                        (end.x, end.y),
                    )
                    .with_stops(&*stops),
                )
            }
        };

        render_context.set_fill_rule(vello_cpu::peniko::Fill::EvenOdd);
        render_context.set_transform(transform);
        render_context.set_paint(paint);
        render_context.fill_path(&path);
        render_context
            .render_to_pixmap(pixels, vello_cpu::RenderMode::default());

        if border_width > 0.0 {
            // Border path is offset by half the border width
            let border_bounds = Rectangle {
                x: quad.bounds.x + border_width / 2.0,
                y: quad.bounds.y + border_width / 2.0,
                width: quad.bounds.width - border_width,
                height: quad.bounds.height - border_width,
            };

            // Make sure the border radius is correct
            let border_radius = <[f32; 4]>::from(quad.border.radius);
            let border_path = rounded_rectangle(border_bounds, border_radius);

            // Stroking a path works well in this case
            render_context.set_paint(vello_cpu::PaintType::Solid(into_color(
                quad.border.color,
            )));
            render_context
                .set_stroke(vello_cpu::kurbo::Stroke::new(border_width.into()));
            render_context.set_transform(transform);
            render_context.fill_path(&border_path);
            render_context
                .render_to_pixmap(pixels, vello_cpu::RenderMode::default());
        }
    }

    pub fn draw_text(
        &mut self,
        text: &Text,
        transformation: Transformation,
        pixels: &mut vello_cpu::Pixmap,
        render_context: &mut vello_cpu::RenderContext,
        clip_bounds: Rectangle,
    ) {
        match text {
            Text::Paragraph {
                paragraph,
                position,
                color,
                clip_bounds: _, // TODO
                transformation: local_transformation,
            } => {
                let transformation = transformation * *local_transformation;

                let physical_bounds =
                    Rectangle::new(*position, paragraph.min_bounds)
                        * transformation;

                if !clip_bounds.intersects(&physical_bounds) {
                    return;
                }
                dbg!(clip_bounds, physical_bounds);

                self.text_pipeline.draw_paragraph(
                    paragraph,
                    *position,
                    *color,
                    pixels,
                    render_context,
                    transformation,
                );
            }
            Text::Editor {
                editor,
                position,
                color,
                clip_bounds: _, // TODO
                transformation: local_transformation,
            } => {
                let transformation = transformation * *local_transformation;

                let physical_bounds =
                    Rectangle::new(*position, editor.bounds) * transformation;

                if !clip_bounds.intersects(&physical_bounds) {
                    return;
                }

                self.text_pipeline.draw_editor(
                    editor,
                    *position,
                    *color,
                    pixels,
                    render_context,
                    transformation,
                );
            }
            Text::Cached {
                content,
                bounds,
                color,
                size,
                line_height,
                font,
                align_x,
                align_y,
                shaping,
                clip_bounds: text_bounds, // TODO
            } => {
                let physical_bounds = *text_bounds * transformation;

                if !clip_bounds.intersects(&physical_bounds) {
                    return;
                }

                self.text_pipeline.draw_cached(
                    content,
                    *bounds,
                    *color,
                    *size,
                    *line_height,
                    *font,
                    *align_x,
                    *align_y,
                    *shaping,
                    pixels,
                    render_context,
                    transformation,
                );
            }
            Text::Raw {
                raw,
                transformation: local_transformation,
            } => {
                let Some(buffer) = raw.buffer.upgrade() else {
                    return;
                };

                let transformation = transformation * *local_transformation;
                let (width, height) = buffer.size();

                let physical_bounds = Rectangle::new(
                    raw.position,
                    Size::new(
                        width.unwrap_or(clip_bounds.width),
                        height.unwrap_or(clip_bounds.height),
                    ),
                ) * transformation;

                if !clip_bounds.intersects(&physical_bounds) {
                    return;
                }

                self.text_pipeline.draw_raw(
                    &buffer,
                    raw.position,
                    raw.color,
                    pixels,
                    render_context,
                    transformation,
                );
            }
        }
    }

    #[cfg(feature = "geometry")]
    pub fn draw_primitive(
        &mut self,
        primitive: &Primitive,
        transformation: Transformation,
        pixels: &mut tiny_skia::PixmapMut<'_>,
        clip_mask: &mut tiny_skia::Mask,
        layer_bounds: Rectangle,
    ) {
        match primitive {
            Primitive::Fill { path, paint, rule } => {
                let physical_bounds = {
                    let bounds = path.bounds();

                    Rectangle {
                        x: bounds.x(),
                        y: bounds.y(),
                        width: bounds.width(),
                        height: bounds.height(),
                    } * transformation
                };

                let Some(clip_bounds) =
                    layer_bounds.intersection(&physical_bounds)
                else {
                    return;
                };

                let clip_mask =
                    (physical_bounds != clip_bounds).then_some(clip_mask as &_);

                pixels.fill_path(
                    path,
                    paint,
                    *rule,
                    into_transform(transformation),
                    clip_mask,
                );
            }
            Primitive::Stroke {
                path,
                paint,
                stroke,
            } => {
                let physical_bounds = {
                    let bounds = path.bounds();

                    Rectangle {
                        x: bounds.x(),
                        y: bounds.y(),
                        width: bounds.width(),
                        height: bounds.height(),
                    } * transformation
                };

                let Some(clip_bounds) =
                    layer_bounds.intersection(&physical_bounds)
                else {
                    return;
                };

                let clip_mask =
                    (physical_bounds != clip_bounds).then_some(clip_mask as &_);

                pixels.stroke_path(
                    path,
                    paint,
                    stroke,
                    into_transform(transformation),
                    clip_mask,
                );
            }
        }
    }

    #[cfg(feature = "image")]
    pub fn draw_image(
        &mut self,
        image: &Image,
        _transformation: Transformation,
        _pixels: &mut tiny_skia::PixmapMut<'_>,
        _clip_mask: &mut tiny_skia::Mask,
        _clip_bounds: Rectangle,
    ) {
        match image {
            #[cfg(feature = "image")]
            Image::Raster(raster, bounds) => {
                let physical_bounds = *bounds * _transformation;

                if !_clip_bounds.intersects(&physical_bounds) {
                    return;
                }

                let clip_mask = (!physical_bounds.is_within(&_clip_bounds))
                    .then_some(_clip_mask as &_);

                let center = physical_bounds.center();
                let radians = f32::from(raster.rotation);

                let transform = into_transform(_transformation).post_rotate_at(
                    radians.to_degrees(),
                    center.x,
                    center.y,
                );

                self.raster_pipeline.draw(
                    &raster.handle,
                    raster.filter_method,
                    *bounds,
                    raster.opacity,
                    _pixels,
                    transform,
                    clip_mask,
                );
            }
            #[cfg(feature = "svg")]
            Image::Vector(svg, bounds) => {
                let physical_bounds = *bounds * _transformation;

                if !_clip_bounds.intersects(&physical_bounds) {
                    return;
                }

                let clip_mask = (!physical_bounds.is_within(&_clip_bounds))
                    .then_some(_clip_mask as &_);

                let center = physical_bounds.center();
                let radians = f32::from(svg.rotation);

                let transform = into_transform(_transformation).post_rotate_at(
                    radians.to_degrees(),
                    center.x,
                    center.y,
                );

                self.vector_pipeline.draw(
                    &svg.handle,
                    svg.color,
                    physical_bounds,
                    svg.opacity,
                    _pixels,
                    transform,
                    clip_mask,
                );
            }
            #[cfg(not(feature = "image"))]
            Image::Raster { .. } => {
                log::warn!(
                    "Unsupported primitive in `iced_tiny_skia`: {image:?}",
                );
            }
            #[cfg(not(feature = "svg"))]
            Image::Vector { .. } => {
                log::warn!(
                    "Unsupported primitive in `iced_tiny_skia`: {image:?}",
                );
            }
        }
    }

    pub fn trim(&mut self) {
        self.text_pipeline.trim_cache();

        #[cfg(feature = "image")]
        self.raster_pipeline.trim_cache();

        #[cfg(feature = "svg")]
        self.vector_pipeline.trim_cache();
    }
}

pub fn into_color(
    color: Color,
) -> vello_cpu::peniko::color::AlphaColor<vello_cpu::peniko::color::Srgb> {
    vello_cpu::peniko::Color::new([color.r, color.g, color.b, color.a])
}

fn into_transform(transformation: Transformation) -> vello_cpu::kurbo::Affine {
    let translation = transformation.translation();

    vello_cpu::kurbo::Affine::new([
        transformation.scale_factor().into(),
        0.0,
        0.0,
        transformation.scale_factor().into(),
        translation.x.into(),
        translation.y.into(),
    ])
}

fn rounded_rectangle(
    bounds: Rectangle,
    border_radius: [f32; 4],
) -> vello_cpu::kurbo::BezPath {
    let [top_left, top_right, bottom_right, bottom_left] = border_radius;

    let rect = vello_cpu::kurbo::Rect::new(
        bounds.x.into(),
        bounds.y.into(),
        bounds.width.into(),
        bounds.height.into(),
    );
    rect.to_rounded_rect((
        top_left.into(),
        top_right.into(),
        bottom_right.into(),
        bottom_left.into(),
    ))
    .to_path(0.1)
}

fn maybe_line_to(path: &mut tiny_skia::PathBuilder, x: f32, y: f32) {
    if path.last_point() != Some(tiny_skia::Point { x, y }) {
        path.line_to(x, y);
    }
}

fn arc_to(
    path: &mut tiny_skia::PathBuilder,
    x_from: f32,
    y_from: f32,
    x_to: f32,
    y_to: f32,
    radius: f32,
) {
    let svg_arc = kurbo::SvgArc {
        from: kurbo::Point::new(f64::from(x_from), f64::from(y_from)),
        to: kurbo::Point::new(f64::from(x_to), f64::from(y_to)),
        radii: kurbo::Vec2::new(f64::from(radius), f64::from(radius)),
        x_rotation: 0.0,
        large_arc: false,
        sweep: true,
    };

    match kurbo::Arc::from_svg_arc(&svg_arc) {
        Some(arc) => {
            arc.to_cubic_beziers(0.1, |p1, p2, p| {
                path.cubic_to(
                    p1.x as f32,
                    p1.y as f32,
                    p2.x as f32,
                    p2.y as f32,
                    p.x as f32,
                    p.y as f32,
                );
            });
        }
        None => {
            path.line_to(x_to, y_to);
        }
    }
}

fn smoothstep(a: f32, b: f32, x: f32) -> f32 {
    let x = ((x - a) / (b - a)).clamp(0.0, 1.0);

    x * x * (3.0 - 2.0 * x)
}

fn rounded_box_sdf(
    to_center: Vector,
    size: tiny_skia::Size,
    radii: &[f32],
) -> f32 {
    let radius = match (to_center.x > 0.0, to_center.y > 0.0) {
        (true, true) => radii[2],
        (true, false) => radii[1],
        (false, true) => radii[3],
        (false, false) => radii[0],
    };

    let x = (to_center.x.abs() - size.width() + radius).max(0.0);
    let y = (to_center.y.abs() - size.height() + radius).max(0.0);

    (x.powf(2.0) + y.powf(2.0)).sqrt() - radius
}

pub fn adjust_clip_mask(
    render_context: &mut vello_cpu::RenderContext,
    bounds: Rectangle,
) {
    render_context.reset();

    let rect = vello_cpu::kurbo::Rect::new(
        bounds.x.into(),
        bounds.y.into(),
        bounds.width.into(),
        bounds.height.into(),
    );

    render_context.set_fill_rule(vello_cpu::peniko::Fill::EvenOdd);
    render_context.fill_rect(&rect);
}
