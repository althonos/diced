use std::ops::Deref;
use std::ops::Range;

use super::Crispr;

#[derive(Debug, Clone, PartialEq)]
pub enum RegionType {
    Spacer,
    Repeat,
}

/// A sequence region.
#[derive(Debug)]
pub struct Region<S> {
    sequence: S,
    start: usize,
    end: usize,
}

impl<S> Region<S> {
    #[inline]
    pub fn new(sequence: S, start: usize, end: usize) -> Self {
        Self {
            sequence,
            start,
            end,
        }
    }

    #[inline]
    pub fn is_empty(&self) -> bool {
        self.start == self.end
    }

    #[inline]
    pub fn len(&self) -> usize {
        self.end - self.start
    }

    #[inline]
    pub fn start(&self) -> usize {
        self.start
    }

    #[inline]
    pub fn end(&self) -> usize {
        self.end
    }
}

impl<S: AsRef<str>> Region<S> {
    #[inline]
    pub fn as_str(&self) -> &str {
        self.as_ref()
    }

    #[inline]
    pub fn as_bytes(&self) -> &[u8] {
        self.as_ref().as_bytes()
    }
}

impl<S: AsRef<str>> AsRef<str> for Region<S> {
    #[inline]
    fn as_ref(&self) -> &str {
        &self.sequence.as_ref()[self.start..self.end]
    }
}

impl<S: AsRef<str>> Deref for Region<S> {
    type Target = str;
    #[inline]
    fn deref(&self) -> &Self::Target {
        self.as_ref()
    }
}

impl<S: AsRef<str>> PartialEq<&str> for Region<S> {
    #[inline]
    fn eq(&self, other: &&str) -> bool {
        self.as_ref() == *other
    }
}

/// An iterator over arbitrary sequence regions.
#[derive(Debug)]
pub struct Regions<'c, S> {
    indices: Range<usize>,
    crispr: &'c Crispr<S>,
    ty: RegionType,
}

impl<'c, S: AsRef<str> + Clone> Regions<'c, S> {
    pub(crate) fn new(crispr: &'c Crispr<S>, ty: RegionType) -> Self {
        Self {
            indices: 0..crispr.len(),
            crispr,
            ty,
        }
    }

    #[inline]
    fn get(&self, i: usize) -> Region<S> {
        match self.ty {
            RegionType::Spacer => self.crispr.spacer(i),
            RegionType::Repeat => self.crispr.repeat(i),
        }
    }
}

impl<'c, S: AsRef<str> + Clone> Iterator for Regions<'c, S> {
    type Item = Region<S>;

    #[inline]
    fn next(&mut self) -> Option<Self::Item> {
        self.indices.next().map(|index| self.get(index))
    }

    #[inline]
    fn nth(&mut self, n: usize) -> Option<Self::Item> {
        self.indices.nth(n).map(|index| self.get(index))
    }
}

impl<'c, S: AsRef<str> + Clone> DoubleEndedIterator for Regions<'c, S> {
    #[inline]
    fn next_back(&mut self) -> Option<Self::Item> {
        self.indices.next_back().map(|index| self.get(index))
    }

    #[inline]
    fn nth_back(&mut self, n: usize) -> Option<Self::Item> {
        self.indices.nth_back(n).map(|index| self.get(index))
    }
}

impl<'c, S: AsRef<str> + Clone> ExactSizeIterator for Regions<'c, S> {
    #[inline]
    fn len(&self) -> usize {
        self.indices.len()
    }
}
