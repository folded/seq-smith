use ndarray::{Array1, Array2};
use numpy::PyReadonlyArray2;
use pyo3::prelude::*;
use pyo3::types::PyBytes;
use pyo3::wrap_pyfunction;
use pyo3_stub_gen::{derive::*, define_stub_info_gatherer};

/// Represents the type of an alignment fragment.
#[gen_stub_pyclass_enum]
#[pyclass(eq)]
#[derive(PartialEq, Eq, Debug, Clone, Copy)]
enum FragType {
    AGap = 0,
    BGap = 1,
    Match = 2,
}

/// Represents a single fragment within a sequence alignment.
///
/// Args:
///     frag_type (FragType): The type of the fragment (e.g., Match, AGap, BGap).
///     sa_start (int): The starting position in sequence A.
///     sb_start (int): The starting position in sequence B.
///     len (int): The length of the fragment.
#[gen_stub_pyclass]
#[pyclass]
#[derive(PartialEq, Eq, Debug, Clone)]
struct AlignFrag {
    #[pyo3(get, set)]
    frag_type: FragType,
    #[pyo3(get)]
    sa_start: i32,
    #[pyo3(get)]
    sb_start: i32,
    #[pyo3(get)]
    len: i32,
}

#[pymethods]
impl AlignFrag {
    #[new]
    fn new(frag_type: FragType, sa_start: i32, sb_start: i32, len: i32) -> Self {
        AlignFrag {
            frag_type,
            sa_start,
            sb_start,
            len,
        }
    }

    fn __repr__(&self) -> String {
        format!(
            "AlignFrag(frag_type={:?}, sa_start={}, sb_start={}, len={})",
            self.frag_type, self.sa_start, self.sb_start, self.len
        )
    }

    fn __eq__(&self, other: &Self) -> bool {
        self.frag_type == other.frag_type
            && self.sa_start == other.sa_start
            && self.sb_start == other.sb_start
            && self.len == other.len
    }
}

/// Represents a complete sequence alignment.
///
/// Args:
///     align_frag (list[AlignFrag]): A list of alignment fragments.
///     frag_count (int): The number of fragments in the alignment.
///     score (int): The total score of the alignment.
#[gen_stub_pyclass]
#[pyclass]
#[derive(Debug, Clone)]
struct Alignment {
    #[pyo3(get)]
    align_frag: Vec<AlignFrag>,
    #[pyo3(get)]
    frag_count: usize,
    #[pyo3(get)]
    score: i32,
}

struct AlignmentParams {
    sa: Vec<u8>,
    sb: Vec<u8>,
    sa_len: usize,
    sb_len: usize,
    score_matrix: Array2<i32>,
    gap_open: i32,
    gap_extend: i32,
}

impl AlignmentParams {
    fn new(
        seqa: Vec<u8>,
        seqb: Vec<u8>,
        score_matrix: Array2<i32>,
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
            sa_len: seqa.len(),
            sb_len: seqb.len(),
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
                curr_score: Array1::uninit(params.sb_len).assume_init(),
                prev_score: Array1::uninit(params.sb_len).assume_init(),
                dir_matrix: Array2::uninit((params.sa_len, params.sb_len)).assume_init(),
                hgap_pos: Array1::uninit(params.sb_len).assume_init(),
                hgap_score: Array1::uninit(params.sb_len).assume_init(),
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
        if score.saturating_add(params.gap_open) >= self.vgap_score.saturating_add(params.gap_extend) {
            self.vgap_score = score.saturating_add(params.gap_open);
            self.vgap_pos = row as i32;
        } else {
            self.vgap_score = self.vgap_score.saturating_add(params.gap_extend);
        }
        if score.saturating_add(params.gap_open) >= self.hgap_score[row].saturating_add(params.gap_extend) {
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
    s_col: usize,
    s_row: usize,
    global_a: bool,
    global_b: bool,
) -> Vec<AlignFrag> {
    let mut result = Vec::new();
    let mut s_col = s_col as i32;
    let mut s_row = s_row as i32;
    let mut d_kind = DirectionKind::Match;

    while s_col >= 0 && s_row >= 0 {
        d_kind = data.dir_matrix[[s_col as usize, s_row as usize]].kind();

        let mut temp = AlignFrag {
            frag_type: FragType::Match,
            sa_start: 0,
            sb_start: 0,
            len: 0,
        };

        match d_kind {
            DirectionKind::Stop => break,
            DirectionKind::GapA(len) => {
                s_row -= len;
                temp.frag_type = FragType::AGap;
                temp.len = len;
            }
            DirectionKind::GapB(len) => {
                s_col -= len;
                temp.frag_type = FragType::BGap;
                temp.len = len;
            }
            DirectionKind::Match => {
                let mut count = 0;
                loop {
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
                temp.frag_type = FragType::Match;
                temp.len = count;
            }
        }
        temp.sa_start = s_col + 1;
        temp.sb_start = s_row + 1;
        result.push(temp);
    }

    if !matches!(d_kind, DirectionKind::Stop) {
        if global_b && s_row >= 0 {
            result.push(AlignFrag {
                frag_type: FragType::AGap,
                sa_start: 0,
                sb_start: 0,
                len: s_row + 1,
            });
        }
        if global_a && s_col >= 0 {
            result.push(AlignFrag {
                frag_type: FragType::BGap,
                sa_start: 0,
                sb_start: 0,
                len: s_col + 1,
            });
        }
    }
    result.reverse();
    result
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

    for row in 0..params.sb_len {
        data.hgap_pos[row] = -1;
        data.hgap_score[row] = params.gap_open;
        data.prev_score[row] = 0;
    }

    for col in 0..params.sa_len {
        data.vgap_pos = -1;
        data.vgap_score = params.gap_open;

        let (score, dir) = data.compute_cell_clipped(0, col, params.match_score(0, col));
        update_max_score(score, 0, col);
        data.write_cell(0, col, score, dir);
        data.update_gaps(0, col, score, &params);

        for row in 1..params.sb_len {
            let match_score = data.prev_score[row - 1].saturating_add(params.match_score(row, col));
            let (score, dir) = data.compute_cell_clipped(row, col, match_score);
            update_max_score(score, row, col);
            data.write_cell(row, col, score, dir);
            data.update_gaps(row, col, score, &params);
        }
        data.swap_scores();
    }

    let align_frag = traceback(&data, max_col, max_row, false, false);

    Ok(Alignment {
        frag_count: align_frag.len(),
        align_frag: align_frag,
        score: max_score,
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
    let params = AlignmentParams::new(seqa.as_bytes().to_vec(), seqb.as_bytes().to_vec(), score_matrix.as_array().into_owned(), gap_open, gap_extend)?;

    py.detach(move || {
        _local_align_core(params)
    })
}

fn _global_align_core(params: AlignmentParams) -> PyResult<Alignment> {
    let mut data = AlignmentData::new(&params);

    for row in 0..params.sb_len {
        let score = params.gap_cost(row as i32 + 1);
        data.prev_score[row] = score;
        data.hgap_pos[row] = -1;
        data.hgap_score[row] = score.saturating_add(params.gap_open);
    }

    for col in 0..params.sa_len {
        data.vgap_pos = -1;
        data.vgap_score = params
            .gap_cost(col as i32 + 1)
            .saturating_add(params.gap_open);

        let match_score = params
            .match_score(0, col)
            .saturating_add(params.gap_cost(col as i32));

        let (score, _) = data.compute_and_write_cell(0, col, match_score);
        data.update_gaps(0, col, score, &params);

        for row in 1..params.sb_len {
            let match_score = data.prev_score[row - 1].saturating_add(params.match_score(row, col));
            let (score, _) = data.compute_and_write_cell(row, col, match_score);
            data.update_gaps(row, col, score, &params);
        }
        data.swap_scores();
    }

    let final_score = data.prev_score[params.sb_len - 1];
    let align_frag = traceback(&data, params.sa_len - 1, params.sb_len - 1, true, true);

    Ok(Alignment {
        frag_count: align_frag.len(),
        align_frag: align_frag,
        score: final_score,
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
    let params = AlignmentParams::new(seqa.as_bytes().to_vec(), seqb.as_bytes().to_vec(), score_matrix.as_array().into_owned(), gap_open, gap_extend)?;

    py.detach(move || {
        _global_align_core(params)
    })
}

fn _local_global_align_core(params: AlignmentParams) -> PyResult<Alignment> {
    let mut data = AlignmentData::new(&params);

    let mut max_score = std::i32::MIN;
    let mut max_row = 0;
    let mut max_col = 0;

    for row in 0..params.sb_len {
        data.prev_score[row] = std::i32::MIN;
        data.hgap_pos[row] = -1;
        data.hgap_score[row] = std::i32::MIN;
    }

    for col in 0..params.sa_len {
        data.vgap_pos = -1;
        data.vgap_score = params.gap_open;

        let (score, _) = data.compute_and_write_cell(0, col, params.match_score(0, col));
        data.update_gaps(0, col, score, &params);

        for row in 1..params.sb_len {
            let match_score = data.prev_score[row - 1].saturating_add(params.match_score(row, col));
            let (score, _) = data.compute_and_write_cell(row, col, match_score);
            data.update_gaps(row, col, score, &params);
        }

        if data.curr_score[params.sb_len - 1] >= max_score {
            max_row = params.sb_len - 1;
            max_col = col;
            max_score = data.curr_score[params.sb_len - 1];
        }
        data.swap_scores();
    }

    let align_frag = traceback(&data, max_col, max_row, false, true);

    Ok(Alignment {
        frag_count: align_frag.len(),
        align_frag: align_frag,
        score: max_score,
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
    let params = AlignmentParams::new(seqa.as_bytes().to_vec(), seqb.as_bytes().to_vec(), score_matrix.as_array().into_owned(), gap_open, gap_extend)?;

    py.detach(move || {
        _local_global_align_core(params)
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

    for row in 0..params.sb_len {
        data.prev_score[row] = 0;
        data.hgap_pos[row] = -1;
        data.hgap_score[row] = params.gap_open;
    }

    for col in 0..params.sa_len {
        data.vgap_pos = -1;
        data.vgap_score = params.gap_open;

        let (score, dir) = data.compute_cell(0, col, params.match_score(0, col));

        data.write_cell(0, col, score, dir);
        data.update_gaps(0, col, score, &params);

        for row in 1..params.sb_len {
            let match_score = data.prev_score[row - 1].saturating_add(params.match_score(row, col));
            let (score, dir) = data.compute_cell(row, col, match_score);

            data.write_cell(row, col, score, dir);
            data.update_gaps(row, col, score, &params);
        }
        update_max_score(data.curr_score[params.sb_len - 1], params.sb_len - 1, col);
        data.swap_scores();
    }
    for row in 0..params.sb_len {
        update_max_score(data.prev_score[row], row, params.sa_len - 1);
    }

    if max_score == std::i32::MIN {
        return Ok(Alignment {
            align_frag: Vec::new(),
            frag_count: 0,
            score: 0,
        });
    }
    let align_frag = traceback(&data, max_col as usize, max_row as usize, false, false);

    Ok(Alignment {
        frag_count: align_frag.len(),
        align_frag: align_frag,
        score: max_score,
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
    let params = AlignmentParams::new(seqa.as_bytes().to_vec(), seqb.as_bytes().to_vec(), score_matrix.as_array().into_owned(), gap_open, gap_extend)?;

    py.detach(move || {
        _overlap_align_core(params)
    })
}

#[pymodule]
fn _seq_smith(_py: Python, m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_wrapped(wrap_pyfunction!(local_align))?;
    m.add_wrapped(wrap_pyfunction!(global_align))?;
    m.add_wrapped(wrap_pyfunction!(local_global_align))?;
    m.add_wrapped(wrap_pyfunction!(overlap_align))?;
    m.add_class::<Alignment>()?;
    m.add_class::<AlignFrag>()?;
    m.add_class::<FragType>()?;
    Ok(())
}

// Define a function to gather stub information
define_stub_info_gatherer!(stub_info);
