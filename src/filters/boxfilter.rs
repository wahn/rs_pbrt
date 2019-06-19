//std
use std::sync::Arc;
// pbrt
use crate::core::filter::Filter;
use crate::core::geometry::{Point2f, Vector2f};
use crate::core::paramset::ParamSet;
use crate::core::pbrt::Float;

// see box.h

#[derive(Debug, Default, Copy, Clone)]
pub struct BoxFilter {
    // inherited from Filter (see filter.h)
    pub radius: Vector2f,
    pub inv_radius: Vector2f,
}

impl BoxFilter {
    pub fn create(ps: &ParamSet) -> Arc<Filter + Sync + Send> {
        let xw: Float = ps.find_one_float("xwidth", 0.5);
        let yw: Float = ps.find_one_float("ywidth", 0.5);
        let box_filter: Arc<Filter + Sync + Send> = Arc::new(BoxFilter {
            radius: Vector2f { x: xw, y: yw },
            inv_radius: Vector2f {
                x: 1.0 / xw,
                y: 1.0 / yw,
            },
        });
        box_filter
    }
}

impl Filter for BoxFilter {
    fn evaluate(&self, _p: Point2f) -> Float {
        1.0
    }
    fn get_radius(&self) -> Vector2f {
        Vector2f {
            x: self.radius.x,
            y: self.radius.y,
        }
    }
}
