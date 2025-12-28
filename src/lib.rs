use numpy::ndarray::{Array1, Array2};
use numpy::PyReadonlyArray2;
use pyo3::prelude::*;
use pyo3::types::PyBytes;
use pyo3::wrap_pyfunction;
use pyo3_stub_gen::{define_stub_info_gatherer, derive::*};
use rayon::prelude::*;
use std::cmp::Ordering;
use std::collections::BinaryHeap;

/// Represents the type of an alignment fragment.
#[gen_stub_pyclass_enum]
#[pyclass(eq)]
#[derive(PartialEq, Eq, Debug, Clone, Copy)]
enum FragmentType {
    AGap = 0,
    BGap = 1,
    Match = 2,
}

/// Represents a single fragment within a sequence alignment.
///
/// Args:
///     fragment_type (FragmentType): The type of the fragment (e.g., Match, AGap, BGap).
///     sa_start (int): The starting position in sequence A.
///     sb_start (int): The starting position in sequence B.
///     len (int): The length of the fragment.
#[gen_stub_pyclass]
#[pyclass]
#[derive(PartialEq, Eq, Debug, Clone)]
struct AlignmentFragment {
    #[pyo3(get, set)]
    fragment_type: FragmentType,
    #[pyo3(get)]
    sa_start: i32,
    #[pyo3(get)]
    sb_start: i32,
    #[pyo3(get)]
    len: i32,
}

#[pymethods]
impl AlignmentFragment {
    #[new]
    fn new(fragment_type: FragmentType, sa_start: i32, sb_start: i32, len: i32) -> Self {
        AlignmentFragment {
            fragment_type,
            sa_start,
            sb_start,
            len,
        }
    }

    fn __repr__(&self) -> String {
        format!(
            "AlignmentFragment(fragment_type={:?}, sa_start={}, sb_start={}, len={})",
            self.fragment_type, self.sa_start, self.sb_start, self.len
        )
    }

    fn __eq__(&self, other: &Self) -> bool {
        self.fragment_type == other.fragment_type
            && self.sa_start == other.sa_start
            && self.sb_start == other.sb_start
            && self.len == other.len
    }
}

/// Represents detailed statistics about an alignment.
#[gen_stub_pyclass]
#[pyclass]
#[derive(Debug, Clone, Copy, Default)]
struct AlignmentStats {
    #[pyo3(get)]
    num_exact_matches: i32,
    #[pyo3(get)]
    num_positive_mismatches: i32,
    #[pyo3(get)]
    num_negative_mismatches: i32,
    #[pyo3(get)]
    num_a_gaps: i32,
    #[pyo3(get)]
    num_b_gaps: i32,
}

#[pymethods]
impl AlignmentStats {
    #[getter]
    fn len(&self) -> i32 {
        self.num_exact_matches
            + self.num_positive_mismatches
            + self.num_negative_mismatches
            + self.num_a_gaps
            + self.num_b_gaps
    }
}

/// Represents a complete sequence alignment.
///
/// Args:
///     fragments (list[AlignmentFragment]): A list of alignment fragments.
///     score (int): The total score of the alignment.
///     stats (AlignmentStats): Statistics about the alignment.
#[gen_stub_pyclass]
#[pyclass]
#[derive(Debug, Clone)]
struct Alignment {
    #[pyo3(get)]
    fragments: Vec<AlignmentFragment>,
    #[pyo3(get)]
    score: i32,
    #[pyo3(get)]
    stats: AlignmentStats,
}

struct AlignmentParams<'a> {
    sa: &'a Vec<u8>,
    sb: &'a Vec<u8>,
    score_matrix: &'a Array2<i32>,
    gap_open: i32,
    gap_extend: i32,
}

impl<'a> AlignmentParams<'a> {
    fn new(
        seqa: &'a Vec<u8>,
        seqb: &'a Vec<u8>,
        score_matrix: &'a Array2<i32>,
        gap_open: i32,
        gap_extend: i32,
    ) -> PyResult<Self> {
        if seqa.is_empty() || seqb.is_empty() {
            return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
                "Input sequences cannot be empty.",
            ));
        }
        if score_matrix.ndim() != 2 {
            return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
                "Score matrix must be 2-dimensional.",
            ));
        }
        let (rows, cols) = score_matrix.dim();
        if rows != cols {
            return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
                "Score matrix must be square.",
            ));
        }
        if gap_open >= 0 {
            return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
                "Gap open penalty must be negative.",
            ));
        }
        if gap_extend >= 0 {
            return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
                "Gap extend penalty must be negative.",
            ));
        }
        Ok(Self {
            sa: seqa,
            sb: seqb,
            score_matrix,
            gap_open,
            gap_extend,
        })
    }

    #[inline(always)]
    fn match_score(&self, row: usize, col: usize) -> i32 {
        self.score_matrix[[self.sa[col] as usize, self.sb[row] as usize]]
    }

    #[inline(always)]
    fn gap_cost(&self, gap_len: i32) -> i32 {
        if gap_len == 0 {
            0
        } else {
            self.gap_open
                .saturating_add(self.gap_extend.saturating_mul(gap_len - 1))
        }
    }
}

struct UngappedAlignmentParams<'a> {
    sa: &'a Vec<u8>,
    sb: &'a Vec<u8>,
    score_matrix: &'a Array2<i32>,
}

impl<'a> UngappedAlignmentParams<'a> {
    fn new(seqa: &'a Vec<u8>, seqb: &'a Vec<u8>, score_matrix: &'a Array2<i32>) -> PyResult<Self> {
        if seqa.is_empty() || seqb.is_empty() {
            return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
                "Input sequences cannot be empty.",
            ));
        }
        if score_matrix.ndim() != 2 {
            return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
                "Score matrix must be 2-dimensional.",
            ));
        }
        let (rows, cols) = score_matrix.dim();
        if rows != cols {
            return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
                "Score matrix must be square.",
            ));
        }
        Ok(Self {
            sa: seqa,
            sb: seqb,
            score_matrix,
        })
    }

    #[inline(always)]
    fn match_score(&self, row: usize, col: usize) -> i32 {
        self.score_matrix[[self.sa[col] as usize, self.sb[row] as usize]]
    }
}

struct AlignmentData {
    curr_score: Array1<i32>,
    prev_score: Array1<i32>,
    dir_matrix: Array2<Direction>,
    hgap_pos: Array1<i32>,
    hgap_score: Array1<i32>,
    vgap_pos: i32,
    vgap_score: i32,
}

impl AlignmentData {
    fn new(params: &AlignmentParams) -> Self {
        unsafe {
            Self {
                curr_score: Array1::uninit(params.sb.len()).assume_init(),
                prev_score: Array1::uninit(params.sb.len()).assume_init(),
                dir_matrix: Array2::uninit((params.sa.len(), params.sb.len())).assume_init(),
                hgap_pos: Array1::uninit(params.sb.len()).assume_init(),
                hgap_score: Array1::uninit(params.sb.len()).assume_init(),
                vgap_pos: -1,
                vgap_score: 0,
            }
        }
    }

    #[allow(dead_code)]
    fn debug_col_scores(&self, col: usize) {
        let sb_len = self.curr_score.len();
        let mut line = format!("[{}]", col);
        for row in 0..sb_len {
            let dir_char = match self.dir_matrix[[col, row]].kind() {
                DirectionKind::Match => "\\",
                DirectionKind::GapA(_) => "-",
                DirectionKind::GapB(_) => "|",
                DirectionKind::Stop => "*",
            };
            let score = self.curr_score[row];
            line.push_str(&format!(" {} {:<3}", dir_char, score));
        }
        eprintln!("{}", line);
    }

    #[inline(always)]
    fn swap_scores(&mut self) {
        std::mem::swap(&mut self.curr_score, &mut self.prev_score);
    }

    #[inline(always)]
    fn compute_cell(&self, row: usize, col: usize, mut score: i32) -> (i32, Direction) {
        let mut dir = Direction::MATCH;
        if score < self.vgap_score {
            score = self.vgap_score;
            dir = Direction::gap_a((row as i32).saturating_sub(self.vgap_pos));
        }
        if score < self.hgap_score[row] {
            score = self.hgap_score[row];
            dir = Direction::gap_b((col as i32).saturating_sub(self.hgap_pos[row]));
        }
        (score, dir)
    }

    #[inline(always)]
    fn compute_cell_clipped(&self, row: usize, col: usize, score: i32) -> (i32, Direction) {
        let (score, dir) = self.compute_cell(row, col, score);
        if score < 0 {
            (0, Direction::STOP)
        } else {
            (score, dir)
        }
    }

    #[inline(always)]
    fn update_gaps(&mut self, row: usize, col: usize, score: i32, params: &AlignmentParams) {
        if score.saturating_add(params.gap_open)
            >= self.vgap_score.saturating_add(params.gap_extend)
        {
            self.vgap_score = score.saturating_add(params.gap_open);
            self.vgap_pos = row as i32;
        } else {
            self.vgap_score = self.vgap_score.saturating_add(params.gap_extend);
        }
        if score.saturating_add(params.gap_open)
            >= self.hgap_score[row].saturating_add(params.gap_extend)
        {
            self.hgap_score[row] = score.saturating_add(params.gap_open);
            self.hgap_pos[row] = col as i32;
        } else {
            self.hgap_score[row] = self.hgap_score[row].saturating_add(params.gap_extend);
        }
    }

    #[inline(always)]
    fn write_cell(&mut self, row: usize, col: usize, score: i32, dir: Direction) {
        self.dir_matrix[[col, row]] = dir;
        self.curr_score[row] = score;
    }

    #[inline(always)]
    fn compute_and_write_cell(
        &mut self,
        row: usize,
        col: usize,
        match_score: i32,
    ) -> (i32, Direction) {
        let (score, dir) = self.compute_cell(row, col, match_score);
        self.write_cell(row, col, score, dir);
        (score, dir)
    }
}

#[derive(Clone, Copy, Debug)]
struct Direction(i32);

enum DirectionKind {
    Match,
    Stop,
    GapA(i32),
    GapB(i32),
}

impl Direction {
    const MATCH: Self = Self(0);
    const STOP: Self = Self(i32::MIN);

    #[inline(always)]
    fn gap_a(len: i32) -> Self {
        debug_assert!(len > 0);
        Self(-len)
    }

    #[inline(always)]
    fn gap_b(len: i32) -> Self {
        debug_assert!(len > 0);
        Self(len)
    }

    #[inline(always)]
    fn kind(&self) -> DirectionKind {
        match self.0 {
            0 => DirectionKind::Match,
            i32::MIN => DirectionKind::Stop,
            val if val < 0 => DirectionKind::GapA(-val),
            val => DirectionKind::GapB(val),
        }
    }
}

fn traceback(
    data: &AlignmentData,
    params: &AlignmentParams,
    s_col: usize,
    s_row: usize,
    global_a: bool,
    global_b: bool,
) -> (Vec<AlignmentFragment>, AlignmentStats) {
    let mut result = Vec::new();
    let mut stats = AlignmentStats::default();

    let mut s_col = s_col as i32;
    let mut s_row = s_row as i32;
    let mut d_kind = DirectionKind::Match;

    while s_col >= 0 && s_row >= 0 {
        d_kind = data.dir_matrix[[s_col as usize, s_row as usize]].kind();

        let mut temp = AlignmentFragment {
            fragment_type: FragmentType::Match,
            sa_start: 0,
            sb_start: 0,
            len: 0,
        };

        match d_kind {
            DirectionKind::Stop => break,
            DirectionKind::GapA(len) => {
                s_row -= len;
                temp.fragment_type = FragmentType::AGap;
                temp.len = len;
                stats.num_a_gaps += len;
            }
            DirectionKind::GapB(len) => {
                s_col -= len;
                temp.fragment_type = FragmentType::BGap;
                temp.len = len;
                stats.num_b_gaps += len;
            }
            DirectionKind::Match => {
                let mut count = 0;
                loop {
                    // Statistics calculation
                    if params.sa[s_col as usize] == params.sb[s_row as usize] {
                        stats.num_exact_matches += 1;
                    } else {
                        let score = params.match_score(s_row as usize, s_col as usize);
                        if score > 0 {
                            stats.num_positive_mismatches += 1;
                        } else {
                            stats.num_negative_mismatches += 1;
                        }
                    }

                    s_col -= 1;
                    s_row -= 1;
                    count += 1;
                    if s_col < 0 || s_row < 0 {
                        break;
                    }
                    if !matches!(
                        data.dir_matrix[[s_col as usize, s_row as usize]].kind(),
                        DirectionKind::Match
                    ) {
                        break;
                    }
                }
                temp.fragment_type = FragmentType::Match;
                temp.len = count;
            }
        }
        temp.sa_start = s_col + 1;
        temp.sb_start = s_row + 1;
        result.push(temp);
    }

    if !matches!(d_kind, DirectionKind::Stop) {
        if global_b && s_row >= 0 {
            result.push(AlignmentFragment {
                fragment_type: FragmentType::AGap,
                sa_start: 0,
                sb_start: 0,
                len: s_row + 1,
            });
            stats.num_a_gaps += s_row + 1;
        }
        if global_a && s_col >= 0 {
            result.push(AlignmentFragment {
                fragment_type: FragmentType::BGap,
                sa_start: 0,
                sb_start: 0,
                len: s_col + 1,
            });
            stats.num_b_gaps += s_col + 1;
        }
    }
    result.reverse();
    (result, stats)
}

fn _local_align_core(params: AlignmentParams) -> PyResult<Alignment> {
    let mut data = AlignmentData::new(&params);

    let mut max_score = 0;
    let mut max_row = 0;
    let mut max_col = 0;

    let mut update_max_score = |score: i32, row: usize, col: usize| {
        if score >= max_score {
            max_row = row;
            max_col = col;
            max_score = score;
        }
    };

    for row in 0..params.sb.len() {
        data.hgap_pos[row] = -1;
        data.hgap_score[row] = params.gap_open;
        data.prev_score[row] = 0;
    }

    for col in 0..params.sa.len() {
        data.vgap_pos = -1;
        data.vgap_score = params.gap_open;

        let (score, dir) = data.compute_cell_clipped(0, col, params.match_score(0, col));
        update_max_score(score, 0, col);
        data.write_cell(0, col, score, dir);
        data.update_gaps(0, col, score, &params);

        for row in 1..params.sb.len() {
            let match_score = data.prev_score[row - 1].saturating_add(params.match_score(row, col));
            let (score, dir) = data.compute_cell_clipped(row, col, match_score);
            update_max_score(score, row, col);
            data.write_cell(row, col, score, dir);
            data.update_gaps(row, col, score, &params);
        }
        data.swap_scores();
    }

    let (fragments, stats) = traceback(&data, &params, max_col, max_row, false, false);

    Ok(Alignment {
        fragments: fragments,
        score: max_score,
        stats,
    })
}

/// Performs a local alignment between two sequences using the Smith-Waterman algorithm.
///
/// Args:
///     seqa (bytes): The first sequence as a byte array.
///     seqb (bytes): The second sequence as a byte array.
///     score_matrix (numpy.ndarray): A 2D numpy array representing the scoring matrix.
///     gap_open (int): The penalty for opening a gap. Must be negative.
///     gap_extend (int): The penalty for extending a gap. Must be negative.
///
/// Raises:
///     ValueError: If any of the following are true:
///         * input sequences are empty
///         * gap penalties are not negative.
///         * score matrix is not 2-dimensional and square.
///
/// Returns:
///     Alignment: An Alignment object containing the score and alignment fragments.
#[gen_stub_pyfunction]
#[pyfunction]
fn local_align<'py>(
    py: Python<'py>,
    seqa: &Bound<'py, PyBytes>,
    seqb: &Bound<'py, PyBytes>,
    score_matrix: PyReadonlyArray2<i32>,
    gap_open: i32,
    gap_extend: i32,
) -> PyResult<Alignment> {
    let seqa = seqa.as_bytes().to_vec();
    let seqb = seqb.as_bytes().to_vec();
    let score_matrix = score_matrix.as_array().into_owned();

    py.detach(move || {
        let params = AlignmentParams::new(&seqa, &seqb, &score_matrix, gap_open, gap_extend)?;
        _local_align_core(params)
    })
}

/// Performs local alignment of one sequence against many sequences in parallel.
#[gen_stub_pyfunction]
#[pyfunction]
#[pyo3(signature = (seqa, seqbs, score_matrix, gap_open, gap_extend, num_threads=None))]
fn local_align_many<'py>(
    py: Python<'py>,
    seqa: &Bound<'py, PyBytes>,
    seqbs: Vec<Bound<'py, PyBytes>>,
    score_matrix: PyReadonlyArray2<i32>,
    gap_open: i32,
    gap_extend: i32,
    num_threads: Option<usize>,
) -> PyResult<Vec<Alignment>> {
    let seqa = seqa.as_bytes().to_vec();
    let seqbs: Vec<Vec<u8>> = seqbs.iter().map(|s| s.as_bytes().to_vec()).collect();
    let score_matrix = score_matrix.as_array().into_owned();

    py.detach(move || {
        _align_many_core(
            seqa,
            seqbs,
            score_matrix,
            gap_open,
            gap_extend,
            num_threads,
            _local_align_core,
        )
    })
}

fn _global_align_core(params: AlignmentParams) -> PyResult<Alignment> {
    let mut data = AlignmentData::new(&params);

    for row in 0..params.sb.len() {
        let score = params.gap_cost(row as i32 + 1);
        data.prev_score[row] = score;
        data.hgap_pos[row] = -1;
        data.hgap_score[row] = score.saturating_add(params.gap_open);
    }

    for col in 0..params.sa.len() {
        data.vgap_pos = -1;
        data.vgap_score = params
            .gap_cost(col as i32 + 1)
            .saturating_add(params.gap_open);

        let match_score = params
            .match_score(0, col)
            .saturating_add(params.gap_cost(col as i32));

        let (score, _) = data.compute_and_write_cell(0, col, match_score);
        data.update_gaps(0, col, score, &params);

        for row in 1..params.sb.len() {
            let match_score = data.prev_score[row - 1].saturating_add(params.match_score(row, col));
            let (score, _) = data.compute_and_write_cell(row, col, match_score);
            data.update_gaps(row, col, score, &params);
        }
        data.swap_scores();
    }

    let final_score = data.prev_score[params.sb.len() - 1];
    let (fragments, stats) = traceback(
        &data,
        &params,
        params.sa.len() - 1,
        params.sb.len() - 1,
        true,
        true,
    );

    Ok(Alignment {
        fragments: fragments,
        score: final_score,
        stats,
    })
}

/// Performs a global alignment between two sequences using the Needleman-Wunsch algorithm.
///
/// Args:
///     seqa (bytes): The first sequence as a byte array.
///     seqb (bytes): The second sequence as a byte array.
///     score_matrix (numpy.ndarray): A 2D numpy array representing the scoring matrix.
///     gap_open (int): The penalty for opening a gap. Must be negative.
///     gap_extend (int): The penalty for extending a gap. Must be negative.
///
/// Raises:
///     ValueError: If any of the following are true:
///         * input sequences are empty
///         * gap penalties are not negative.
///         * score matrix is not 2-dimensional and square.
///
/// Returns:
///     Alignment: An Alignment object containing the score and alignment fragments.
#[gen_stub_pyfunction]
#[pyfunction]
fn global_align<'py>(
    py: Python<'py>,
    seqa: &Bound<'py, PyBytes>,
    seqb: &Bound<'py, PyBytes>,
    score_matrix: PyReadonlyArray2<i32>,
    gap_open: i32,
    gap_extend: i32,
) -> PyResult<Alignment> {
    let seqa = seqa.as_bytes().to_vec();
    let seqb = seqb.as_bytes().to_vec();
    let score_matrix = score_matrix.as_array().into_owned();

    py.detach(move || {
        let params = AlignmentParams::new(&seqa, &seqb, &score_matrix, gap_open, gap_extend)?;
        _global_align_core(params)
    })
}

/// Performs global alignment of one sequence against many sequences in parallel.
#[gen_stub_pyfunction]
#[pyfunction]
#[pyo3(signature = (seqa, seqbs, score_matrix, gap_open, gap_extend, num_threads=None))]
fn global_align_many<'py>(
    py: Python<'py>,
    seqa: &Bound<'py, PyBytes>,
    seqbs: Vec<Bound<'py, PyBytes>>,
    score_matrix: PyReadonlyArray2<i32>,
    gap_open: i32,
    gap_extend: i32,
    num_threads: Option<usize>,
) -> PyResult<Vec<Alignment>> {
    let seqa = seqa.as_bytes().to_vec();
    let seqbs: Vec<Vec<u8>> = seqbs.iter().map(|s| s.as_bytes().to_vec()).collect();
    let score_matrix = score_matrix.as_array().into_owned();

    py.detach(move || {
        _align_many_core(
            seqa,
            seqbs,
            score_matrix,
            gap_open,
            gap_extend,
            num_threads,
            _global_align_core,
        )
    })
}

fn _local_global_align_core(params: AlignmentParams) -> PyResult<Alignment> {
    let mut data = AlignmentData::new(&params);

    let mut max_score = std::i32::MIN;
    let mut max_row = 0;
    let mut max_col = 0;

    for row in 0..params.sb.len() {
        let score = params.gap_cost(row as i32 + 1);
        data.prev_score[row] = score;
        data.hgap_pos[row] = -1;
        data.hgap_score[row] = score.saturating_add(params.gap_open);
    }

    for col in 0..params.sa.len() {
        data.vgap_pos = -1;
        data.vgap_score = params.gap_open;

        let (score, _) = data.compute_and_write_cell(0, col, params.match_score(0, col));
        data.update_gaps(0, col, score, &params);

        for row in 1..params.sb.len() {
            let match_score = data.prev_score[row - 1].saturating_add(params.match_score(row, col));
            let (score, _) = data.compute_and_write_cell(row, col, match_score);
            data.update_gaps(row, col, score, &params);
        }

        if data.curr_score[params.sb.len() - 1] >= max_score {
            max_row = params.sb.len() - 1;
            max_col = col;
            max_score = data.curr_score[params.sb.len() - 1];
        }
        data.swap_scores();
    }

    let (fragments, stats) = traceback(&data, &params, max_col, max_row, false, true);

    Ok(Alignment {
        fragments: fragments,
        score: max_score,
        stats,
    })
}

/// Performs a local-global alignment. This alignment finds the best local alignment of `seqa`
/// within `seqb`, but `seqb` must be aligned globally.
///
/// Args:
///     seqa (bytes): The first sequence as a byte array.
///     seqb (bytes): The second sequence as a byte array.
///     score_matrix (numpy.ndarray): A 2D numpy array representing the scoring matrix.
///     gap_open (int): The penalty for opening a gap. Must be negative.
///     gap_extend (int): The penalty for extending a gap. Must be negative.
///
/// Raises:
///     ValueError: If any of the following are true:
///         * input sequences are empty
///         * gap penalties are not negative.
///         * score matrix is not 2-dimensional and square.
///
/// Returns:
///     Alignment: An Alignment object containing the score and alignment fragments.
#[gen_stub_pyfunction]
#[pyfunction]
fn local_global_align<'py>(
    py: Python<'py>,
    seqa: &Bound<'py, PyBytes>,
    seqb: &Bound<'py, PyBytes>,
    score_matrix: PyReadonlyArray2<i32>,
    gap_open: i32,
    gap_extend: i32,
) -> PyResult<Alignment> {
    let seqa = seqa.as_bytes().to_vec();
    let seqb = seqb.as_bytes().to_vec();
    let score_matrix = score_matrix.as_array().into_owned();

    py.detach(move || {
        let params = AlignmentParams::new(&seqa, &seqb, &score_matrix, gap_open, gap_extend)?;
        _local_global_align_core(params)
    })
}

/// Performs local-global alignment of one sequence against many sequences in parallel.
#[gen_stub_pyfunction]
#[pyfunction]
#[pyo3(signature = (seqa, seqbs, score_matrix, gap_open, gap_extend, num_threads=None))]
fn local_global_align_many<'py>(
    py: Python<'py>,
    seqa: &Bound<'py, PyBytes>,
    seqbs: Vec<Bound<'py, PyBytes>>,
    score_matrix: PyReadonlyArray2<i32>,
    gap_open: i32,
    gap_extend: i32,
    num_threads: Option<usize>,
) -> PyResult<Vec<Alignment>> {
    let seqa = seqa.as_bytes().to_vec();
    let seqbs: Vec<Vec<u8>> = seqbs.iter().map(|s| s.as_bytes().to_vec()).collect();
    let score_matrix = score_matrix.as_array().into_owned();

    py.detach(move || {
        _align_many_core(
            seqa,
            seqbs,
            score_matrix,
            gap_open,
            gap_extend,
            num_threads,
            _local_global_align_core,
        )
    })
}

fn _overlap_align_core(params: AlignmentParams) -> PyResult<Alignment> {
    // An overlap alignment must start on the bottom or right edge of the DP matrix.
    // Gaps at the start are not penalized.

    let mut data = AlignmentData::new(&params);

    let mut max_score = std::i32::MIN;
    let mut max_row = -1;
    let mut max_col = -1;

    let mut update_max_score = |score: i32, row: usize, col: usize| {
        if score >= max_score {
            max_row = row as i32;
            max_col = col as i32;
            max_score = score;
        }
    };

    for row in 0..params.sb.len() {
        data.prev_score[row] = 0;
        data.hgap_pos[row] = -1;
        data.hgap_score[row] = params.gap_open;
    }

    for col in 0..params.sa.len() {
        data.vgap_pos = -1;
        data.vgap_score = params.gap_open;

        let (score, dir) = data.compute_cell(0, col, params.match_score(0, col));

        data.write_cell(0, col, score, dir);
        data.update_gaps(0, col, score, &params);

        for row in 1..params.sb.len() {
            let match_score = data.prev_score[row - 1].saturating_add(params.match_score(row, col));
            let (score, dir) = data.compute_cell(row, col, match_score);

            data.write_cell(row, col, score, dir);
            data.update_gaps(row, col, score, &params);
        }
        update_max_score(
            data.curr_score[params.sb.len() - 1],
            params.sb.len() - 1,
            col,
        );
        data.swap_scores();
    }
    for row in 0..params.sb.len() {
        update_max_score(data.prev_score[row], row, params.sa.len() - 1);
    }

    if max_score == std::i32::MIN {
        return Ok(Alignment {
            fragments: Vec::new(),
            score: 0,
            stats: AlignmentStats::default(),
        });
    }
    let (fragments, stats) = traceback(
        &data,
        &params,
        max_col as usize,
        max_row as usize,
        false,
        false,
    );

    Ok(Alignment {
        fragments: fragments,
        score: max_score,
        stats,
    })
}

/// Performs an overlap alignment between two sequences.
///
/// This alignment type does not penalize gaps at the start or end of either sequence,
/// making it suitable for finding overlaps between sequences.
///
/// Args:
///     seqa (bytes): The first sequence as a byte array.
///     seqb (bytes): The second sequence as a byte array.
///     score_matrix (numpy.ndarray): A 2D numpy array representing the scoring matrix.
///     gap_open (int): The penalty for opening a gap. Must be negative.
///     gap_extend (int): The penalty for extending a gap. Must be negative.
///
/// Raises:
///     ValueError: If any of the following are true:
///         * input sequences are empty
///         * gap penalties are not negative.
///         * score matrix is not 2-dimensional and square.
///
/// Returns:
///     Alignment: An Alignment object containing the score and alignment fragments.
#[gen_stub_pyfunction]
#[pyfunction]
fn overlap_align<'py>(
    py: Python<'py>,
    seqa: &Bound<'py, PyBytes>,
    seqb: &Bound<'py, PyBytes>,
    score_matrix: PyReadonlyArray2<i32>,
    gap_open: i32,
    gap_extend: i32,
) -> PyResult<Alignment> {
    let seqa = seqa.as_bytes().to_vec();
    let seqb = seqb.as_bytes().to_vec();
    let score_matrix = score_matrix.as_array().into_owned();

    py.detach(move || {
        let params = AlignmentParams::new(&seqa, &seqb, &score_matrix, gap_open, gap_extend)?;
        _overlap_align_core(params)
    })
}

/// Performs overlap alignment of one sequence against many sequences in parallel.
#[gen_stub_pyfunction]
#[pyfunction]
#[pyo3(signature = (seqa, seqbs, score_matrix, gap_open, gap_extend, num_threads=None))]
fn overlap_align_many<'py>(
    py: Python<'py>,
    seqa: &Bound<'py, PyBytes>,
    seqbs: Vec<Bound<'py, PyBytes>>,
    score_matrix: PyReadonlyArray2<i32>,
    gap_open: i32,
    gap_extend: i32,
    num_threads: Option<usize>,
) -> PyResult<Vec<Alignment>> {
    let seqa = seqa.as_bytes().to_vec();
    let seqbs: Vec<Vec<u8>> = seqbs.iter().map(|s| s.as_bytes().to_vec()).collect();
    let score_matrix = score_matrix.as_array().into_owned();

    py.detach(move || {
        _align_many_core(
            seqa,
            seqbs,
            score_matrix,
            gap_open,
            gap_extend,
            num_threads,
            _overlap_align_core,
        )
    })
}

fn _align_many_core<F>(
    seqa: Vec<u8>,
    seqbs: Vec<Vec<u8>>,
    score_matrix: Array2<i32>,
    gap_open: i32,
    gap_extend: i32,
    num_threads: Option<usize>,
    align_func: F,
) -> PyResult<Vec<Alignment>>
where
    F: Fn(AlignmentParams) -> PyResult<Alignment> + Sync + Send,
{
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads.unwrap_or(0)) // 0 tells rayon to use a default number of threads
        .build()
        .map_err(|e| {
            PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!(
                "Failed to create thread pool: {}",
                e
            ))
        })?;

    pool.install(|| {
        seqbs
            .into_par_iter()
            .map(|seqb| {
                let params =
                    AlignmentParams::new(&seqa, &seqb, &score_matrix, gap_open, gap_extend)?;
                align_func(params)
            })
            .collect()
    })
}

#[derive(Eq, PartialEq)]
struct Candidate {
    score: i32,
    sa_start: usize,
    sb_start: usize,
    len: usize,
}

impl Ord for Candidate {
    fn cmp(&self, other: &Self) -> Ordering {
        self.score
            .cmp(&other.score)
            .then_with(|| self.sa_start.cmp(&other.sa_start))
            .then_with(|| self.sb_start.cmp(&other.sb_start))
    }
}

impl PartialOrd for Candidate {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

fn _top_k_ungapped_local_align_core(
    params: UngappedAlignmentParams,
    k: usize,
    filter_overlap_a: bool,
    filter_overlap_b: bool,
) -> PyResult<Vec<Alignment>> {
    let sa_len = params.sa.len();
    let sb_len = params.sb.len();

    let mut candidates: BinaryHeap<Candidate> = BinaryHeap::new();

    let mut add_candidate = |score: i32, sa_start: usize, sb_start: usize, len: usize| {
        if score > 0 {
            candidates.push(Candidate {
                score,
                sa_start,
                sb_start,
                len,
            });
        }
    };

    let mut process_diagonal = |start_row: usize, start_col: usize, max_len: usize| {
        let mut curr_score = 0;
        let mut segment_start_idx = 0; // index along diagonal where current positive segment started
        let mut peak_score = 0;
        let mut peak_idx = 0; // index along diagonal where peak occurred

        for i in 0..max_len {
            let row = start_row + i;
            let col = start_col + i;
            let val = params.match_score(row, col);

            if curr_score == 0 && val <= 0 {
                continue;
            }
            if curr_score == 0 {
                segment_start_idx = i;
            }

            curr_score += val;

            if curr_score <= 0 {
                add_candidate(
                    peak_score,
                    start_col + segment_start_idx,
                    start_row + segment_start_idx,
                    peak_idx - segment_start_idx + 1,
                );
                curr_score = 0;
                peak_score = 0;
            } else {
                if curr_score > peak_score {
                    peak_score = curr_score;
                    peak_idx = i;
                }
            }
        }
        add_candidate(
            peak_score,
            start_col + segment_start_idx,
            start_row + segment_start_idx,
            peak_idx - segment_start_idx + 1,
        );
    };

    // Diagonals starting at first row (row=0, col=0..sa_len)
    for start_col in 0..sa_len {
        let max_len = std::cmp::min(sa_len - start_col, sb_len);
        process_diagonal(0, start_col, max_len);
    }

    // Diagonals starting at first column (row=1..sb_len, col=0)
    for start_row in 1..sb_len {
        let max_len = std::cmp::min(sa_len, sb_len - start_row);
        process_diagonal(start_row, 0, max_len);
    }

    // Select top k non-overlapping
    let mut alignments: Vec<Alignment> = Vec::with_capacity(k);

    while alignments.len() < k {
        if let Some(candidate) = candidates.pop() {
            // Check overlap
            let sa_end = candidate.sa_start + candidate.len;
            let sb_end = candidate.sb_start + candidate.len;

            let overlap = alignments.iter().any(|prev| {
                let p_sa_start = prev.fragments[0].sa_start as usize;
                let p_sb_start = prev.fragments[0].sb_start as usize;
                let p_sa_end = p_sa_start + prev.fragments[0].len as usize;
                let p_sb_end = p_sb_start + prev.fragments[0].len as usize;

                // Overlap in A?
                let overlaps_a =
                    filter_overlap_a && candidate.sa_start < p_sa_end && sa_end > p_sa_start;
                // Overlap in B?
                let overlaps_b =
                    filter_overlap_b && candidate.sb_start < p_sb_end && sb_end > p_sb_start;

                overlaps_a || overlaps_b
            });

            if !overlap {
                // Construct Alignment
                let mut stats = AlignmentStats::default();
                for i in 0..candidate.len {
                    let r = candidate.sb_start + i;
                    let c = candidate.sa_start + i;
                    if params.sb[r] == params.sa[c] {
                        stats.num_exact_matches += 1;
                    } else if params.match_score(r, c) > 0 {
                        stats.num_positive_mismatches += 1;
                    } else {
                        stats.num_negative_mismatches += 1;
                    }
                }

                alignments.push(Alignment {
                    fragments: vec![AlignmentFragment {
                        fragment_type: FragmentType::Match,
                        sa_start: candidate.sa_start as i32,
                        sb_start: candidate.sb_start as i32,
                        len: candidate.len as i32,
                    }],
                    score: candidate.score,
                    stats: stats,
                });
            }
        } else {
            break;
        }
    }

    Ok(alignments)
}

/// Finds the top-k non-overlapping ungapped local alignments (HSPs).
///
/// Args:
///     seqa (bytes): The first sequence.
///     seqb (bytes): The second sequence.
///     score_matrix (numpy.ndarray): Scorin matrix.
///     k (int): Number of alignments to return.
///
/// Returns:
///     list[Alignment]: List of top-k non-overlapping alignments.
#[gen_stub_pyfunction]
#[pyfunction]
#[pyo3(signature = (seqa, seqb, score_matrix, k, filter_overlap_a=true, filter_overlap_b=true))]
fn top_k_ungapped_local_align<'py>(
    py: Python<'py>,
    seqa: &Bound<'py, PyBytes>,
    seqb: &Bound<'py, PyBytes>,
    score_matrix: PyReadonlyArray2<i32>,
    k: usize,
    filter_overlap_a: bool,
    filter_overlap_b: bool,
) -> PyResult<Vec<Alignment>> {
    let seqa = seqa.as_bytes().to_vec();
    let seqb = seqb.as_bytes().to_vec();
    let score_matrix = score_matrix.as_array().into_owned();

    py.detach(move || {
        _top_k_ungapped_local_align_core(
            UngappedAlignmentParams::new(&seqa, &seqb, &score_matrix)?,
            k,
            filter_overlap_a,
            filter_overlap_b,
        )
    })
}

/// Finds the top-k non-overlapping ungapped local alignments (HSPs) against many sequences in parallel.
///
/// Args:
///     seqa (bytes): The query sequence.
///     seqbs (list[bytes]): List of target sequences.
///     score_matrix (numpy.ndarray): Scoring matrix.
///     k (int): Number of alignments to return per target sequence.
///     num_threads (int, optional): Number of threads to use. Defaults to all available.
///
/// Returns:
///     list[list[Alignment]]: List of alignment lists.
#[gen_stub_pyfunction]
#[pyfunction]
#[pyo3(signature = (seqa, seqbs, score_matrix, k, num_threads=None, filter_overlap_a=true, filter_overlap_b=true))]
fn top_k_ungapped_local_align_many<'py>(
    py: Python<'py>,
    seqa: &Bound<'py, PyBytes>,
    seqbs: Vec<Bound<'py, PyBytes>>,
    score_matrix: PyReadonlyArray2<i32>,
    k: usize,
    num_threads: Option<usize>,
    filter_overlap_a: bool,
    filter_overlap_b: bool,
) -> PyResult<Vec<Vec<Alignment>>> {
    let seqa = seqa.as_bytes().to_vec();
    let seqbs: Vec<Vec<u8>> = seqbs.iter().map(|s| s.as_bytes().to_vec()).collect();
    let score_matrix = score_matrix.as_array().into_owned();

    py.detach(move || {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(num_threads.unwrap_or(0))
            .build()
            .map_err(|e| {
                PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!(
                    "Failed to create thread pool: {}",
                    e
                ))
            })?;

        pool.install(|| {
            seqbs
                .into_par_iter()
                .map(|seqb| {
                    _top_k_ungapped_local_align_core(
                        UngappedAlignmentParams::new(&seqa, &seqb, &score_matrix)?,
                        k,
                        filter_overlap_a,
                        filter_overlap_b,
                    )
                })
                .collect()
        })
    })
}

#[pymodule]
fn _seq_smith(_py: Python, m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_wrapped(wrap_pyfunction!(local_align))?;
    m.add_wrapped(wrap_pyfunction!(global_align))?;
    m.add_wrapped(wrap_pyfunction!(local_global_align))?;
    m.add_wrapped(wrap_pyfunction!(overlap_align))?;
    m.add_wrapped(wrap_pyfunction!(local_align_many))?;
    m.add_wrapped(wrap_pyfunction!(global_align_many))?;
    m.add_wrapped(wrap_pyfunction!(local_global_align_many))?;
    m.add_wrapped(wrap_pyfunction!(overlap_align_many))?;
    m.add_wrapped(wrap_pyfunction!(top_k_ungapped_local_align))?;
    m.add_wrapped(wrap_pyfunction!(top_k_ungapped_local_align_many))?;
    m.add_class::<Alignment>()?;
    m.add_class::<AlignmentFragment>()?;
    m.add_class::<FragmentType>()?;
    Ok(())
}

// Define a function to gather stub information
define_stub_info_gatherer!(stub_info);
