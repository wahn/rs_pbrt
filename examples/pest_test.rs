#[macro_use]
extern crate pest;

use pest::prelude::*;
use std::collections::LinkedList;
    
#[derive(Debug, PartialEq)]
pub enum Node {
    Sentence(LinkedList<Node>),
    Word(LinkedList<Node>),
    Letter(char),
}

impl_rdp! {
    grammar! {
        sentence = _{ word ~ ([" "] ~ word)* }
        word     =  { letter* }
        letter   =  { ['a'..'z'] }
    }

    process! {
        main(&self) -> Node {
            (list: _sentence()) => {
                Node::Sentence(list)
            }
        }

        _sentence(&self) -> LinkedList<Node> {
            (_: word, head: _word(), mut tail: _sentence()) => {
                tail.push_front(Node::Word(head));

                tail
            },
            () => {
                LinkedList::new()
            }
        }

        _word(&self) -> LinkedList<Node> {
            (&head: letter, mut tail: _word()) => {
                tail.push_front(Node::Letter(head.chars().next().unwrap()));

                tail
            },
            () => {
                LinkedList::new()
            }
        }
    }
}

fn main() {
    let mut parser = Rdp::new(StringInput::new("abc def"));

    assert!(parser.sentence());
    println!("{:?}", parser.main());
}
