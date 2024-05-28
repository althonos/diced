#[derive(Debug, Clone, Copy, PartialEq)]
enum Flank {
    Left,
    Right,
}

#[derive(Clone)]
pub struct CrisprFinder {
    min_repeat_count: usize,
    min_repeat_length: usize,
    max_repeat_length: usize,
    min_spacer_length: usize,
    max_spacer_length: usize,
    search_window_length: usize,
}

impl CrisprFinder {
    pub fn find_crispr<S: AsRef<str>>(&self, sequence: S) -> CrisprIterator<S> {
        let sequence_length = sequence.as_ref().len();
        let parameters = self.clone();
        CrisprIterator {
            sequence,
            sequence_length,
            parameters,
            j: 0,
        }
    }
}

impl Default for CrisprFinder {
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

pub struct CrisprIterator<S> {
    parameters: CrisprFinder,
    sequence: S,
    sequence_length: usize,
    j: usize,
}

impl<S: AsRef<str>> CrisprIterator<S> {
    const THRESHOLD: f32 = 0.75;
    const SPACER_TO_SPACER_MAX_SIMILARITY: f32 = 0.62;
    const SPACER_TO_SPACER_LENGTH_DIFF: usize = 12;
    const SPACER_TO_REPEAT_LENGTH_DIFF: usize = 30;

    fn _similarity(s1: &str, s2: &str) -> f32 {
        let max_len = s1.len().max(s2.len());
        1.0 - ((strsim::levenshtein(s1, s2) as f32) / (max_len as f32))
    }

    fn _hamming(s1: &str, s2: &str) -> usize {
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

        let mut last_repeat_index = crispr.indices[num_repeats - 1];
        let mut second_to_last_repeat_index = crispr.indices[num_repeats - 2];
        let mut repeat_spacing = last_repeat_index - second_to_last_repeat_index;

        loop {
            let candidate_repeat_index = last_repeat_index + repeat_spacing;
            let mut begin_search = candidate_repeat_index - scan_range;
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

            if let Some(k) = seq[begin_search..end_search].find(pattern) {
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
        let mut char_counts = vec![0; u8::MAX as usize];

        let mut right_extension_length = self.parameters.search_window_length;
        let max_right_extension_length =
            shortest_repeat_spacing - self.parameters.min_spacer_length;

        while right_extension_length <= max_right_extension_length {
            if last_repeat_start_index + right_extension_length >= sequence_len {
                if crispr.indices.len() - 1 > self.parameters.min_repeat_count {
                    crispr.indices.pop().unwrap();
                    last_repeat_start_index = *crispr.indices.last().unwrap();
                } else {
                    break;
                }
            }
            for k in 0..crispr.indices.len() {
                let current_repeat_start_index = crispr.indices[k];
                let current_repeat = &seq[current_repeat_start_index
                    ..current_repeat_start_index + right_extension_length];
                let first_char = *current_repeat.as_bytes().last().unwrap();
                char_counts[first_char as usize] += 1;
            }
            if char_counts
                .iter()
                .any(|&n| ((n as f32) / (crispr.indices.len() as f32)) >= Self::THRESHOLD)
            {
                right_extension_length += 1;
                char_counts.fill(0);
            } else {
                break;
            }
        }
        right_extension_length -= 1;
        char_counts.fill(0);

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
                char_counts[first_char as usize] += 1;
            }
            if char_counts.iter().any(|&n| {
                ((n as f32) / (crispr.indices.len() as f32)) >= CrisprIterator::<S>::THRESHOLD
            }) {
                left_extension_length += 1;
                char_counts.fill(0);
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
                if Self::_similarity(current_spacer, next_spacer)
                    > Self::SPACER_TO_SPACER_MAX_SIMILARITY
                {
                    return false;
                }
                if Self::_similarity(current_repeat, current_spacer)
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
                    (crispr.repeat_spacing(0, 1) + crispr.repeat_spacing(1, 2)) / 2
                } else {
                    crispr.repeat_spacing(0, 1)
                };
                candidate_repeat_index = first_repeat_index.checked_sub(repeat_spacing)?;
            }
            Flank::Right => {
                repeat_string = crispr.repeat(num_repeats - 1);
                let repeat_spacing = if num_repeats >= 3 {
                    (crispr.repeat_spacing(num_repeats - 2, num_repeats - 1)
                        + crispr.repeat_spacing(num_repeats - 3, num_repeats - 2))
                        / 2
                } else {
                    crispr.repeat_spacing(num_repeats - 2, num_repeats - 1)
                };
                candidate_repeat_index = last_repeat_index + repeat_spacing;
            }
        };

        let mut begin = candidate_repeat_index - scan_range;
        let mut end = candidate_repeat_index + scan_range;

        let scan_left_max_end = first_repeat_index
            .checked_sub(repeat_length)?
            .checked_sub(self.parameters.min_spacer_length)?;
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

        // if begin < 0 {
        //     return None;
        // }
        if begin + repeat_length > sequence_len {
            return None;
        }
        if end + repeat_length > sequence_len {
            end = sequence_len - repeat_length;
        }
        if begin >= end {
            return None;
        }

        let mut array = Vec::new();
        for i in begin..=end {
            let candidate_repeat_string = &seq[i..i + repeat_length];
            array.push(Self::_hamming(repeat_string, candidate_repeat_string));
        }

        let new_candidate_repeat_index =
            begin + (0..array.len()).min_by_key(|&i| array[i]).unwrap();
        let new_candidate_repeat_string =
            &seq[new_candidate_repeat_index..new_candidate_repeat_index + repeat_length];
        if Self::_similarity(repeat_string, new_candidate_repeat_string) >= confidence {
            Some(new_candidate_repeat_index)
        } else {
            None
        }
    }

    fn _trim(&self, crispr: &mut Crispr<S>) {
        let num_repeats = crispr.indices.len();

        let mut char_counts = vec![0u32; u8::MAX as usize];
        while num_repeats > self.parameters.min_repeat_length + 1 {
            for k in 0..num_repeats {
                let repeat = crispr.repeat(k).as_bytes();
                let last_char = repeat[repeat.len() - 1];
                char_counts[last_char as usize] += 1;
            }
            if char_counts
                .iter()
                .all(|&n| ((n as f32) / (crispr.indices.len() as f32)) < Self::THRESHOLD)
            {
                crispr.repeat_length -= 1;
                char_counts.fill(0);
            } else {
                break;
            }
        }

        char_counts.fill(0);
        while num_repeats > self.parameters.min_repeat_length + 1 {
            for k in 0..num_repeats {
                let repeat = crispr.repeat(k).as_bytes();
                let first_char = repeat[0];
                char_counts[first_char as usize] += 1;
            }
            if char_counts
                .iter()
                .all(|&n| ((n as f32) / (crispr.indices.len() as f32)) < Self::THRESHOLD)
            {
                for index in crispr.indices.iter_mut() {
                    *index += 1;
                }
                crispr.repeat_length -= 1;
                char_counts.fill(0);
            } else {
                break;
            }
        }
    }
}

impl<S: AsRef<str> + Clone> Iterator for CrisprIterator<S> {
    type Item = Crispr<S>;
    fn next(&mut self) -> Option<Self::Item> {
        let seq = self.sequence.as_ref();

        let skips = self
            .parameters
            .min_repeat_length
            .checked_sub(2 * self.parameters.search_window_length - 1)
            .unwrap_or(1);

        let search_end = seq.len()
            - self.parameters.max_repeat_length
            - self.parameters.max_spacer_length
            - self.parameters.search_window_length;

        while self.j <= search_end {
            let mut candidate_crispr = Crispr::new(self.sequence.clone());

            let begin_search =
                self.j + self.parameters.min_repeat_length + self.parameters.min_spacer_length;
            let mut end_search = self.j
                + self.parameters.max_repeat_length
                + self.parameters.max_spacer_length
                + self.parameters.search_window_length;
            if end_search > self.sequence_length {
                end_search = self.sequence_length;
            } else if end_search < begin_search {
                end_search = begin_search;
            }

            let pattern = &seq[self.j..self.j + self.parameters.search_window_length];
            if let Some(k) = seq[begin_search..end_search].find(pattern) {
                candidate_crispr.indices.push(self.j);
                candidate_crispr.indices.push(begin_search + k);
                self._scan_right(&mut candidate_crispr, pattern, 24);
            }

            if candidate_crispr.indices.len() >= self.parameters.min_repeat_count {
                self._get_actual_repeat_length(&mut candidate_crispr);
                let actual_repeat_length = candidate_crispr.repeat_length;
                if (actual_repeat_length >= self.parameters.min_repeat_length)
                    && (actual_repeat_length <= self.parameters.max_repeat_length)
                    && self._has_non_repeating_spacers(&candidate_crispr)
                    && self._has_similarly_sized_spacers(&candidate_crispr)
                {
                    self._check_flank(&mut candidate_crispr, Flank::Left, 30, 0.7);
                    self._check_flank(&mut candidate_crispr, Flank::Right, 30, 0.7);
                    self._trim(&mut candidate_crispr);
                    self.j = candidate_crispr.end();
                    return Some(candidate_crispr);
                }
            }

            self.j += skips;
        }

        None
    }
}

#[derive(Debug)]
pub struct Crispr<S> {
    sequence: S,
    indices: Vec<usize>,
    repeat_length: usize,
}

impl<S> Crispr<S> {
    pub fn start(&self) -> usize {
        self.indices.first().cloned().unwrap_or(0)
    }

    pub fn end(&self) -> usize {
        self.indices.last().cloned().unwrap_or(0) + self.repeat_length
    }
}

impl<S: AsRef<str>> Crispr<S> {
    pub fn new(sequence: S) -> Self {
        Self {
            sequence,
            indices: Vec::new(),
            repeat_length: 0,
        }
    }

    pub fn region(&self) -> &str {
        &self.sequence.as_ref()[self.start()..self.end()]
    }

    pub fn repeat(&self, index: usize) -> &str {
        let s = self.sequence.as_ref();
        let start = self.indices[index];
        let end = start + self.repeat_length;
        &s[start..end]
    }

    pub fn spacer(&self, index: usize) -> &str {
        let s = self.sequence.as_ref();
        let current_end = self.indices[index] + self.repeat_length;
        let next_start = self.indices[index + 1];
        // let spacer_start = current_end + 1;
        let spacer_start = current_end;
        let spacer_end = next_start;
        &s[spacer_start..spacer_end]
    }

    pub fn repeat_spacing(&self, pos1: usize, pos2: usize) -> usize {
        self.indices[pos1].abs_diff(self.indices[pos2])
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn find_crisprs() {
        const SEQ: &str = concat!(
            "TTTTACAATCTGCGTTTTAACTCCACACGGTACATTAGAAACCATCTGCAACATATT",
            "CAAGTTCAGCTTCAAAACCTTGTTTTAACTCCACACGGTACATTAGAAACTTCGTCA",
            "AGCTTTACCTCAAAAGTCCTCTCAAACCTGTTTTAACTCCACACGGTACATTAGAAA",
            "CAATAATCAACAACTCTTTGATTTTGTGAAATGGAAGAAGTTTTAACTCCACACGGT",
            "ACATTAGAAACAGAACTCTCAGAAGAACCGAGAGCTTTTTCTATTAACGTTTTAACT",
            "CCACACGGTACATTAGAAACCCTGCGTGCCTGTGTCTAAAAAATA",
        );

        let it = CrisprFinder::default().find_crispr(SEQ);
        let crisprs = it.collect::<Vec<_>>();
        assert_eq!(crisprs.len(), 1);

        assert_eq!(crisprs[0].indices.len(), 5);
        assert_eq!(crisprs[0].repeat(0), "GTTTTAACTCCACACGGTACATTAGAAAC");
        assert_eq!(crisprs[0].start(), 13);
        assert_eq!(crisprs[0].end(), 305);

        let region = crisprs[0].region();
        assert!(region.starts_with(crisprs[0].repeat(0)),);
        assert!(region.ends_with(crisprs[0].repeat(4)),);
    }
}
