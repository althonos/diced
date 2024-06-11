#![doc = include_str!("../README.md")]

#[cfg(feature = "memchr")]
extern crate memchr;

mod region;

pub use self::region::Region;
pub use self::region::Regions;

use self::region::RegionType;
use std::ops::Deref;

#[derive(Debug, Clone, Copy, PartialEq)]
enum Flank {
    Left,
    Right,
}

#[derive(Debug, Default)]
struct DnaCount {
    a: usize,
    t: usize,
    c: usize,
    g: usize,
}

impl DnaCount {
    fn new() -> Self {
        Self::default()
    }

    #[inline]
    fn clear(&mut self) {
        self.a = 0;
        self.c = 0;
        self.t = 0;
        self.g = 0;
    }

    #[inline]
    fn count(&mut self, c: char) {
        match c {
            'a' | 'A' => self.a += 1,
            'c' | 'C' => self.c += 1,
            'g' | 'G' => self.g += 1,
            't' | 'T' => self.t += 1,
            _ => (),
        }
    }

    #[inline]
    fn max(&self) -> usize {
        self.a.max(self.c).max(self.t).max(self.g)
    }
}

#[derive(Debug)]
struct Sequence<S> {
    sequence: S,
    mask: Vec<Region<S>>,
}

impl<S: AsRef<str> + Clone> Sequence<S> {
    const MASK_SIZE: usize = 100;
    pub fn new(sequence: S) -> Self {
        let s = sequence.as_ref();
        let mut mask = Vec::new();

        let mut i = 0;
        let mut j;
        let mut n;

        while i < s.len() {
            n = 0;
            j = i + 1;

            while j < s.len() && s.as_bytes()[i] == s.as_bytes()[j] {
                n += 1;
                j += 1;
            }

            if n >= Self::MASK_SIZE {
                mask.push(Region::new(sequence.clone(), i, j));
            }

            i = j;
        }

        // add an empty mask at the end to facilitate some scanning code later
        mask.push(Region::new(sequence.clone(), s.len(), s.len()));

        Self { sequence, mask }
    }

    pub fn len(&self) -> usize {
        self.as_ref().len()
    }

    fn _is_masked(&self, index: &mut usize, begin: usize, end: usize) -> bool {
        while begin > self.mask[*index].end() {
            *index += 1;
        }
        begin <= self.mask[*index].end() && self.mask[*index].start() <= end
    }
}

impl<S> Deref for Sequence<S> {
    type Target = S;
    fn deref(&self) -> &Self::Target {
        &self.sequence
    }
}

/// A builder type to parameterize a [`Scanner`].
#[derive(Clone)]
pub struct ScannerBuilder {
    min_repeat_count: usize,
    min_repeat_length: usize,
    max_repeat_length: usize,
    min_spacer_length: usize,
    max_spacer_length: usize,
    search_window_length: usize,
}

impl ScannerBuilder {
    /// Create a new scanner builder with default parameters.
    pub fn new() -> Self {
        Self::default()
    }

    /// Scan the provided sequence for CRISPR regions iteratively.
    ///
    /// The sequence can be provided as any string view that also implements
    /// [`Clone`]. This allows several smart pointers to be passed (`&str`,
    /// `Rc<str>`, `Arc<str>`, etc.). The actual `sequence` object however
    /// will be cloned into the result [`Crispr`], so make sure it implements
    /// a cheap [`Clone`], and avoid passing a [`String`].
    ///
    pub fn scan<S: AsRef<str> + Clone>(&self, sequence: S) -> Scanner<S> {
        let mut scanner = Scanner::new(sequence);
        self.clone_into(&mut scanner.parameters);
        scanner
    }

    /// Set the minimum repeat number for CRISPR detection.
    pub fn min_repeat_count(&mut self, min_repeat_count: usize) -> &mut Self {
        self.min_repeat_count = min_repeat_count;
        self
    }

    /// Set the minimum length for each CRISPR repeat.
    pub fn min_repeat_length(&mut self, min_repeat_length: usize) -> &mut Self {
        self.min_repeat_length = min_repeat_length;
        self
    }

    /// Set the maximum length for each CRISPR repeat.
    pub fn max_repeat_length(&mut self, max_repeat_length: usize) -> &mut Self {
        self.max_repeat_length = max_repeat_length;
        self
    }

    /// Set the minimum length for each CRISPR spacer.
    pub fn min_spacer_length(&mut self, min_spacer_length: usize) -> &mut Self {
        self.min_spacer_length = min_spacer_length;
        self
    }

    /// Set the maximum length for each CRISPR spacer.
    pub fn max_spacer_length(&mut self, max_spacer_length: usize) -> &mut Self {
        self.max_spacer_length = max_spacer_length;
        self
    }
}

impl Default for ScannerBuilder {
    fn default() -> Self {
        Self {
            min_repeat_count: 3,
            min_repeat_length: 23,
            max_repeat_length: 47,
            min_spacer_length: 26,
            max_spacer_length: 50,
            search_window_length: 8,
        }
    }
}

/// A scanner for identifying CRISPR regions in a nucleotide sequence.
pub struct Scanner<S> {
    parameters: ScannerBuilder,
    sequence: Sequence<S>,
    sequence_length: usize,
    mask_index: usize,
    j: usize,
}

impl<S: AsRef<str> + Clone> Scanner<S> {
    const THRESHOLD: f32 = 0.75;
    const SPACER_TO_SPACER_MAX_SIMILARITY: f32 = 0.62;
    const SPACER_TO_SPACER_LENGTH_DIFF: usize = 12;
    const SPACER_TO_REPEAT_LENGTH_DIFF: usize = 30;

    #[inline]
    pub fn new(sequence: S) -> Self {
        let seq = Sequence::new(sequence);
        Self {
            parameters: ScannerBuilder::default(),
            sequence_length: seq.len(),
            j: 0,
            sequence: seq,
            mask_index: 0,
        }
    }

    #[inline]
    pub fn sequence(&self) -> &S {
        &self.sequence
    }

    fn _similarity<S1: AsRef<str>, S2: AsRef<str>>(s1: S1, s2: S2) -> f32 {
        let s1 = s1.as_ref();
        let s2 = s2.as_ref();
        let max_len = s1.len().max(s2.len());
        1.0 - ((strsim::levenshtein(s1, s2) as f32) / (max_len as f32))
    }

    fn _hamming<S1: AsRef<str>, S2: AsRef<str>>(s1: S1, s2: S2) -> usize {
        let s1 = s1.as_ref();
        let s2 = s2.as_ref();
        if s1.len() == s2.len() {
            strsim::hamming(s1, s2).unwrap()
        } else {
            let l = s1.len().min(s2.len());
            let d = s1.len().abs_diff(s2.len());
            strsim::hamming(&s1[..l], &s2[..l]).unwrap() + d
        }
    }

    fn _scan_right(&self, crispr: &mut Crispr<S>, pattern: &str, scan_range: usize) {
        let seq = crispr.sequence.as_ref();

        let num_repeats = crispr.indices.len();
        let pattern_len = pattern.len();
        let sequence_len = seq.len();

        #[cfg(feature = "memchr")]
        let finder = memchr::memmem::Finder::new(pattern);

        let mut last_repeat_index = crispr.indices[num_repeats - 1];
        let mut second_to_last_repeat_index = crispr.indices[num_repeats - 2];
        let mut repeat_spacing = last_repeat_index - second_to_last_repeat_index;

        loop {
            let candidate_repeat_index = last_repeat_index + repeat_spacing;
            let mut begin_search = candidate_repeat_index.saturating_sub(scan_range);
            let mut end_search = candidate_repeat_index + pattern_len + scan_range;

            let scan_right_min_begin =
                last_repeat_index + pattern_len + self.parameters.min_spacer_length;
            if begin_search < scan_right_min_begin {
                begin_search = scan_right_min_begin;
            }
            if end_search > sequence_len {
                end_search = sequence_len;
            }
            if begin_search >= end_search {
                break;
            }

            let subseq = &seq[begin_search..end_search];

            #[cfg(feature = "memchr")]
            let pos = finder.find(subseq.as_bytes());
            #[cfg(not(feature = "memchr"))]
            let pos = subseq.find(pattern);

            if let Some(k) = pos {
                crispr.indices.push(begin_search + k);
                second_to_last_repeat_index = last_repeat_index;
                last_repeat_index = begin_search + k;
                repeat_spacing = last_repeat_index - second_to_last_repeat_index;
                if repeat_spacing < self.parameters.min_spacer_length + pattern_len {
                    break;
                }
            } else {
                break;
            }
        }
    }

    fn _get_actual_repeat_length(&self, crispr: &mut Crispr<S>) {
        let mut first_repeat_start_index = *crispr.indices.first().unwrap();
        let mut last_repeat_start_index = *crispr.indices.last().unwrap();
        let shortest_repeat_spacing = crispr
            .indices
            .iter()
            .zip(&crispr.indices[1..])
            .map(|(i, j)| j - i)
            .min()
            .unwrap();

        let seq = crispr.sequence.as_ref();
        let sequence_len = seq.len();
        let mut char_counts = DnaCount::new();

        let mut right_extension_length = self.parameters.search_window_length;
        let max_right_extension_length =
            shortest_repeat_spacing - self.parameters.min_spacer_length;

        while right_extension_length <= max_right_extension_length {
            if last_repeat_start_index + right_extension_length >= sequence_len {
                if crispr.indices.len() > self.parameters.min_repeat_count + 1 {
                    crispr.indices.pop().unwrap();
                    last_repeat_start_index = *crispr.indices.last().unwrap();
                } else {
                    break;
                }
            }
            for k in 0..crispr.indices.len() {
                let current_repeat_start_index = crispr.indices[k];
                let last_char =
                    seq.as_bytes()[current_repeat_start_index + right_extension_length - 1];
                char_counts.count(last_char as char);
            }
            if ((char_counts.max() as f32) / (crispr.indices.len() as f32)) >= Self::THRESHOLD {
                right_extension_length += 1;
                char_counts.clear();
            } else {
                break;
            }
        }
        right_extension_length -= 1;
        char_counts.clear();

        let mut left_extension_length = 0;
        let max_left_extension_length =
            shortest_repeat_spacing - self.parameters.min_spacer_length - right_extension_length;
        while left_extension_length <= max_left_extension_length {
            if first_repeat_start_index < left_extension_length {
                if crispr.indices.len() > self.parameters.min_repeat_count + 1 {
                    crispr.indices.remove(0); // FIXME: use VecDeque?
                    first_repeat_start_index = *crispr.indices.first().unwrap();
                } else {
                    break;
                }
            }
            for k in 0..crispr.indices.len() {
                let current_repeat_start_index = crispr.indices[k];
                let first_char = seq.as_bytes()[current_repeat_start_index - left_extension_length];
                char_counts.count(first_char as char)
            }
            if (char_counts.max() as f32) / (crispr.indices.len() as f32) >= Self::THRESHOLD {
                left_extension_length += 1;
                char_counts.clear();
            } else {
                break;
            }
        }
        left_extension_length -= 1;

        for index in crispr.indices.iter_mut() {
            *index -= left_extension_length;
        }

        crispr.repeat_length = right_extension_length + left_extension_length;
    }

    fn _has_non_repeating_spacers(&self, crispr: &Crispr<S>) -> bool {
        let first_repeat = crispr.repeat(0);
        let first_spacer = crispr.spacer(0);

        if crispr.indices.len() >= 3 {
            let mut i = 0;
            while i + 2 < crispr.indices.len() {
                if i == 4 {
                    return true;
                }
                let current_spacer = crispr.spacer(i);
                let next_spacer = crispr.spacer(i + 1);
                let current_repeat = crispr.repeat(i);
                if Self::_similarity(&current_spacer, &next_spacer)
                    > Self::SPACER_TO_SPACER_MAX_SIMILARITY
                {
                    return false;
                }
                if Self::_similarity(&current_repeat, &current_spacer)
                    > Self::SPACER_TO_SPACER_MAX_SIMILARITY
                {
                    return false;
                }
                i += 1;
            }
            Self::_similarity(crispr.repeat(i), crispr.spacer(i))
                <= Self::SPACER_TO_SPACER_MAX_SIMILARITY
        } else if crispr.indices.len() == 2 {
            if first_spacer.is_empty() {
                false
            } else {
                Self::_similarity(first_spacer, first_repeat)
                    < Self::SPACER_TO_SPACER_MAX_SIMILARITY
            }
        } else {
            false
        }
    }

    fn _has_similarly_sized_spacers(&self, crispr: &Crispr<S>) -> bool {
        let initial_spacer_length = crispr.spacer(0).len();
        let repeat_length = crispr.repeat_length;
        for i in 0..crispr.indices.len() - 1 {
            let current_spacer_length = crispr.spacer(i).len();
            if current_spacer_length.abs_diff(initial_spacer_length)
                > Self::SPACER_TO_SPACER_LENGTH_DIFF
            {
                return false;
            }
            if current_spacer_length.abs_diff(repeat_length) > Self::SPACER_TO_REPEAT_LENGTH_DIFF {
                return false;
            }
        }
        true
    }

    fn _check_flank(
        &self,
        crispr: &mut Crispr<S>,
        flank: Flank,
        scan_range: usize,
        confidence: f32,
    ) {
        while let Some(pos) = self._scan(crispr, flank, scan_range, confidence) {
            match flank {
                Flank::Left => crispr.indices.insert(0, pos),
                Flank::Right => crispr.indices.push(pos),
            }
        }
    }

    fn _scan(
        &self,
        crispr: &Crispr<S>,
        flank: Flank,
        scan_range: usize,
        confidence: f32,
    ) -> Option<usize> {
        let repeat_length = crispr.repeat_length;
        let num_repeats = crispr.indices.len();
        let seq = crispr.sequence.as_ref();
        let sequence_len = seq.len();

        let first_repeat_index = *crispr.indices.first().unwrap();
        let last_repeat_index = *crispr.indices.last().unwrap();

        let repeat_string;
        let candidate_repeat_index;

        match flank {
            Flank::Left => {
                repeat_string = crispr.repeat(0);
                let repeat_spacing = if num_repeats >= 3 {
                    (crispr.repeat_spacing(0) + crispr.repeat_spacing(1)) / 2
                } else {
                    crispr.repeat_spacing(0)
                };
                candidate_repeat_index = first_repeat_index.saturating_sub(repeat_spacing);
            }
            Flank::Right => {
                repeat_string = crispr.repeat(num_repeats - 1);
                let repeat_spacing = if num_repeats >= 3 {
                    (crispr.repeat_spacing(num_repeats - 2)
                        + crispr.repeat_spacing(num_repeats - 3))
                        / 2
                } else {
                    crispr.repeat_spacing(num_repeats - 2)
                };
                candidate_repeat_index = last_repeat_index + repeat_spacing;
            }
        };

        if candidate_repeat_index < scan_range {
            return None;
        }

        let mut begin = candidate_repeat_index - scan_range;
        let mut end = candidate_repeat_index + scan_range;

        let scan_left_max_end = first_repeat_index
            .saturating_sub(repeat_length)
            .saturating_sub(self.parameters.min_spacer_length);
        let scan_right_min_begin =
            last_repeat_index + repeat_length + self.parameters.min_spacer_length;

        match flank {
            Flank::Left => {
                if end > scan_left_max_end {
                    end = scan_left_max_end;
                }
            }
            Flank::Right => {
                if begin < scan_right_min_begin {
                    begin = scan_right_min_begin;
                }
            }
        }

        if begin + repeat_length > sequence_len {
            return None;
        }
        if end + repeat_length > sequence_len {
            end = sequence_len.saturating_sub(repeat_length);
        }
        if begin >= end {
            return None;
        }

        let new_candidate_repeat_index = (begin..=end)
            .min_by_key(|&i| Self::_hamming(&repeat_string, &seq[i..i + repeat_length]))
            .unwrap();
        let new_candidate_repeat_string =
            &seq[new_candidate_repeat_index..new_candidate_repeat_index + repeat_length];
        if Self::_similarity(&repeat_string, new_candidate_repeat_string) >= confidence {
            Some(new_candidate_repeat_index)
        } else {
            None
        }
    }

    fn _trim(&self, crispr: &mut Crispr<S>) {
        let num_repeats = crispr.indices.len();

        let mut char_counts = DnaCount::new();
        while crispr.repeat_length > self.parameters.min_repeat_length {
            for k in 0..num_repeats {
                let repeat = crispr.repeat(k);
                let last_char = repeat.as_bytes().last().unwrap();
                char_counts.count(*last_char as char);
            }
            if (char_counts.max() as f32) / (num_repeats as f32) < Self::THRESHOLD {
                crispr.repeat_length -= 1;
                char_counts.clear();
            } else {
                break;
            }
        }

        char_counts.clear();
        while crispr.repeat_length > self.parameters.min_repeat_length {
            for k in 0..num_repeats {
                let repeat = crispr.repeat(k);
                let first_char = repeat.chars().next().unwrap();
                char_counts.count(first_char);
            }
            if (char_counts.max() as f32) / (num_repeats as f32) < Self::THRESHOLD {
                for index in crispr.indices.iter_mut() {
                    *index += 1;
                }
                crispr.repeat_length -= 1;
                char_counts.clear();
            } else {
                break;
            }
        }
    }
}

impl<S: AsRef<str> + Clone> Iterator for Scanner<S> {
    type Item = Crispr<S>;
    fn next(&mut self) -> Option<Self::Item> {
        let seq = self.sequence.as_ref();

        let skips = self
            .parameters
            .min_repeat_length
            .checked_sub(2 * self.parameters.search_window_length - 1)
            .unwrap_or(1);

        let search_end = seq
            .len()
            .saturating_sub(self.parameters.max_repeat_length)
            .saturating_sub(self.parameters.max_spacer_length)
            .saturating_sub(self.parameters.search_window_length);

        while self.j < search_end {
            let mut begin_search =
                self.j + self.parameters.min_repeat_length + self.parameters.min_spacer_length;
            let mut end_search = self.j
                + self.parameters.max_repeat_length
                + self.parameters.max_spacer_length
                + self.parameters.search_window_length;
            if begin_search > self.sequence_length {
                begin_search = self.sequence_length;
            }
            if end_search > self.sequence_length {
                end_search = self.sequence_length;
            }
            if end_search < begin_search {
                end_search = begin_search;
            }

            if self
                .sequence
                ._is_masked(&mut self.mask_index, begin_search, end_search)
            {
                self.j = self.sequence.mask[self.mask_index].end();
                if self.j >= search_end {
                    return None;
                }
                continue;
            }

            let pattern_start = self.j;
            let pattern_end = (self.j + self.parameters.search_window_length).min(seq.len());

            let pattern = &seq[pattern_start..pattern_end];
            let subseq = &seq[begin_search..end_search];

            #[cfg(feature = "memchr")]
            let pos = memchr::memmem::find(subseq.as_bytes(), pattern.as_bytes());
            #[cfg(not(feature = "memchr"))]
            let pos = subseq.find(pattern);

            let mut candidate_crispr = Crispr::new(self.sequence.clone());
            if let Some(k) = pos {
                candidate_crispr.indices.push(self.j);
                candidate_crispr.indices.push(begin_search + k);
                self._scan_right(&mut candidate_crispr, pattern, 24);
            }

            if candidate_crispr.indices.len() >= self.parameters.min_repeat_count {
                self._get_actual_repeat_length(&mut candidate_crispr);
                let actual_repeat_length = candidate_crispr.repeat_length;

                if actual_repeat_length >= self.parameters.min_repeat_length {
                    if actual_repeat_length <= self.parameters.max_repeat_length {
                        if self._has_non_repeating_spacers(&candidate_crispr) {
                            if self._has_similarly_sized_spacers(&candidate_crispr) {
                                self._check_flank(&mut candidate_crispr, Flank::Left, 30, 0.7);
                                self._check_flank(&mut candidate_crispr, Flank::Right, 30, 0.7);
                                self._trim(&mut candidate_crispr);
                                self.j = candidate_crispr.end();
                                return Some(candidate_crispr);
                            }
                        }
                    }
                }
            }

            self.j += skips;
        }

        None
    }
}

/// A CRISPR repeat region in a nucleotide sequence.
#[derive(Debug, Clone)]
pub struct Crispr<S> {
    sequence: S,
    indices: Vec<usize>,
    repeat_length: usize,
}

impl<S> Crispr<S> {
    /// Get the number of repeats in the CRISPR region.
    #[inline]
    pub fn len(&self) -> usize {
        self.indices.len()
    }

    /// Get the start index of the CRISPR region (zero-based).
    ///
    /// This is returned as a zero-based, inclusive index, which can be
    /// used for slicing.
    #[inline]
    pub fn start(&self) -> usize {
        self.indices.first().cloned().unwrap_or(0)
    }

    /// Get the end index of the CRISPR region (zero-based, exclusive).
    #[inline]
    pub fn end(&self) -> usize {
        self.indices.last().cloned().unwrap_or(0) + self.repeat_length
    }
}

impl<S: AsRef<str>> Crispr<S> {
    /// Create a new crispr region for the given sequence.
    #[inline]
    fn new(sequence: S) -> Self {
        Self {
            sequence,
            indices: Vec::new(),
            repeat_length: 0,
        }
    }
}

impl<S: AsRef<str> + Clone> Crispr<S> {
    /// Get the complete CRISPR region as a [`Region`].
    #[inline]
    pub fn to_region(&self) -> Region<S> {
        Region::new(self.sequence.clone(), self.start(), self.end())
    }

    /// Get the sequence of the `k`-th repeat in the CRISPR region.
    ///
    /// # Panic
    /// Panics if `k >= self.len()`.
    pub fn repeat(&self, index: usize) -> Region<S> {
        let start = self.indices[index];
        let end = start + self.repeat_length;
        Region::new(self.sequence.clone(), start, end)
    }

    /// Get an iterator over the repeats of the CRISPR region.
    #[inline]
    pub fn repeats(&self) -> Regions<S> {
        Regions::new(self, RegionType::Repeat)
    }

    /// Get the sequence of the `k`-th repeat in the CRISPR region.
    ///
    /// # Panic
    /// Panics if `k + 1 >= self.len()`.
    pub fn spacer(&self, index: usize) -> Region<S> {
        let current_end = self.indices[index] + self.repeat_length;
        let next_start = self.indices[index + 1];
        let spacer_start = current_end;
        let spacer_end = next_start;
        Region::new(self.sequence.clone(), spacer_start, spacer_end)
    }

    /// Get an iterator over the spacers of the CRISPR region.
    #[inline]
    pub fn spacers(&self) -> Regions<S> {
        Regions::new(self, RegionType::Spacer)
    }

    /// Compute the spacing the `k`-th and `k+1`-th repeats.
    ///
    /// # Panic
    /// Panics if `k + 1 >= self.len()`.
    #[inline]
    fn repeat_spacing(&self, index: usize) -> usize {
        self.indices[index + 1] - self.indices[index]
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    use std::rc::Rc;

    const SEQ: &str = concat!(
        "TTTTACAATCTGCGTTTTAACTCCACACGGTACATTAGAAACCATCTGCAACATATT",
        "CAAGTTCAGCTTCAAAACCTTGTTTTAACTCCACACGGTACATTAGAAACTTCGTCA",
        "AGCTTTACCTCAAAAGTCCTCTCAAACCTGTTTTAACTCCACACGGTACATTAGAAA",
        "CAATAATCAACAACTCTTTGATTTTGTGAAATGGAAGAAGTTTTAACTCCACACGGT",
        "ACATTAGAAACAGAACTCTCAGAAGAACCGAGAGCTTTTTCTATTAACGTTTTAACT",
        "CCACACGGTACATTAGAAACCCTGCGTGCCTGTGTCTAAAAAATA",
    );

    #[test]
    fn scan_str() {
        let it = ScannerBuilder::default().scan(SEQ);
        let crisprs = it.collect::<Vec<_>>();
        assert_eq!(crisprs.len(), 1);

        assert_eq!(crisprs[0].indices.len(), 5);
        assert_eq!(crisprs[0].repeat(0), "GTTTTAACTCCACACGGTACATTAGAAAC");
        assert_eq!(crisprs[0].start(), 13);
        assert_eq!(crisprs[0].end(), 305);

        let region = crisprs[0].to_region();
        assert!(region.starts_with(crisprs[0].repeat(0).as_ref()),);
        assert!(region.ends_with(crisprs[0].repeat(4).as_ref()),);
    }

    #[test]
    fn scan_rc() {
        let it = ScannerBuilder::default().scan(Rc::from(SEQ));
        let crisprs = it.collect::<Vec<_>>();
        assert_eq!(crisprs.len(), 1);

        assert_eq!(crisprs[0].indices.len(), 5);
        assert_eq!(crisprs[0].repeat(0), "GTTTTAACTCCACACGGTACATTAGAAAC");
        assert_eq!(crisprs[0].start(), 13);
        assert_eq!(crisprs[0].end(), 305);

        let region = crisprs[0].to_region();
        assert!(region.starts_with(crisprs[0].repeat(0).as_ref()),);
        assert!(region.ends_with(crisprs[0].repeat(4).as_ref()),);
    }

    #[test]
    fn scan_empty() {
        let it = ScannerBuilder::default().scan("");
        let crisprs = it.collect::<Vec<_>>();
        assert_eq!(crisprs.len(), 0);
    }

    #[test]
    fn scan_max_under_min() {
        let it = ScannerBuilder::default()
            .min_repeat_length(40)
            .max_repeat_length(10)
            .scan(SEQ);
        let crisprs = it.collect::<Vec<_>>();
        assert_eq!(crisprs.len(), 0);

        let it = ScannerBuilder::default()
            .min_spacer_length(40)
            .max_spacer_length(10)
            .scan(SEQ);
        let crisprs = it.collect::<Vec<_>>();
        assert_eq!(crisprs.len(), 0);
    }
}
