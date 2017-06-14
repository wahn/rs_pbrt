#[macro_use]
extern crate pest;

use pest::prelude::*;

impl_rdp! {
    grammar! {
        expression = _{
            { ["("] ~ expression ~ [")"] | number }
            addition       = { plus  | minus } // precedence 0 is addition
            multiplication = { times | slash } // precedence 1 is multiplication
        }
        number = @{ ["-"]? ~ (["0"] | ['1'..'9'] ~ ['0'..'9']*) }
        plus   =  { ["+"] }
        minus  =  { ["-"] }
        times  =  { ["*"] }
        slash  =  { ["/"] }

        whitespace = _{ [" "] }
    }

    process! {
        compute(&self) -> i32 {
            (&number: number) => number.parse::<i32>().unwrap(),
            (_: addition, left: compute(), sign, right: compute()) => {
                match sign.rule {
                    Rule::plus  => left + right,
                    Rule::minus => left - right,
                    _ => unreachable!()
                }
            },
            (_: multiplication, left: compute(), sign, right: compute()) => {
                match sign.rule {
                    Rule::times => left * right,
                    Rule::slash => left / right,
                    _ => unreachable!()
                }
            }
        }
    }
}

fn main() {
    let mut parser = Rdp::new(StringInput::new("(3 + (9 + 3 * 4 + (3 + 1) / 2 - 4)) * 2"));

    assert!(parser.expression());
    assert_eq!(parser.compute(), 44);
}
