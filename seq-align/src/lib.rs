use pyo3::prelude::*;
use pyo3::wrap_pyfunction;
use numpy::{PyReadonlyArray2};

#[pyclass]
#[derive(PartialEq, Eq, Debug, Clone, Copy)]
enum FragType {
    AGap = 0,
    BGap = 1,
    Match = 2,
}

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
        self.frag_type == other.frag_type &&
        self.sa_start == other.sa_start &&
        self.sb_start == other.sb_start &&
        self.len == other.len
    }
}

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

    fn gap_a(len: i32) -> Self {
        debug_assert!(len > 0);
        Self(-len)
    }

    fn gap_b(len: i32) -> Self {
        debug_assert!(len > 0);
        Self(len)
    }

    fn kind(&self) -> DirectionKind {
        match self.0 {
            0 => DirectionKind::Match,
            i32::MIN => DirectionKind::Stop,
            val if val < 0 => DirectionKind::GapA(-val),
            val => DirectionKind::GapB(val),
        }
    }
}


fn traceback(dir_matrix: &[Direction], sb_len: usize, s_col: usize, s_row: usize) -> Vec<AlignFrag> {
    let mut result = Vec::new();
    let mut s_col = s_col as i32;
    let mut s_row = s_row as i32;

    while s_col >= 0 && s_row >= 0 {
        let d_kind = dir_matrix[(s_col as usize) * sb_len + (s_row as usize)].kind();

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
                    if !matches!(dir_matrix[(s_col as usize) * sb_len + (s_row as usize)].kind(), DirectionKind::Match) {
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
    result.reverse();
    result
}

fn global_traceback(
    dir_matrix: &[Direction],
    sb_len: usize,
    s_col: usize,
    s_row: usize,
) -> Vec<AlignFrag> {
    let mut result = Vec::new();
    let mut s_col = s_col as i32;
    let mut s_row = s_row as i32;

    while s_col >= 0 && s_row >= 0 {
        let d_kind = dir_matrix[(s_col as usize) * sb_len + (s_row as usize)].kind();

        let mut temp = AlignFrag {
            frag_type: FragType::Match,
            sa_start: 0,
            sb_start: 0,
            len: 0,
        };

        match d_kind {
            DirectionKind::Stop => break, // Should not happen in global
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
                    if !matches!(dir_matrix[(s_col as usize) * sb_len + (s_row as usize)].kind(), DirectionKind::Match) {
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

    if s_col >= 0 || s_row >= 0 {
        let mut temp = AlignFrag {
            frag_type: FragType::Match,
            sa_start: 0,
            sb_start: 0,
            len: 0,
        };
        temp.sa_start = 0;
        temp.sb_start = 0;
        if s_col >= 0 {
            temp.frag_type = FragType::BGap;
            temp.len = s_col + 1;
        } else {
            temp.frag_type = FragType::AGap;
            temp.len = s_row + 1;
        }
        result.push(temp);
    }
    result.reverse();
    result
}


#[pyfunction]
fn local_align(
    seqa: &[u8],
    seqb: &[u8],
    score_matrix: PyReadonlyArray2<i32>,
    gap_open: i32,
    gap_extend: i32,
) -> PyResult<Alignment> {
    let sa = seqa;
    let sb = seqb;
    let score_matrix = score_matrix.as_array();

    let sa_len = sa.len();
    let sb_len = sb.len();

    if sa_len == 0 || sb_len == 0 {
        return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
            "Input sequences cannot be empty.",
        ));
    }

    let mut curr_score = vec![0; sb_len];
    let mut prev_score = vec![0; sb_len];
    let mut dir_matrix = vec![Direction::STOP; sa_len * sb_len];
    let mut hgap_pos = vec![0; sb_len];
    let mut hgap_score = vec![0; sb_len];

    let mut max_score = 0;
    let mut max_row = 0;
    let mut max_col = 0;

    for row in 0..sb_len {
        hgap_pos[row] = -1;
        hgap_score[row] = gap_open;
        prev_score[row] = 0;
    }

    for col in 0..sa_len {
        let mut vgap_pos = -1;
        let mut vgap_score = gap_open;

        let mut score = score_matrix[[sa[col] as usize, sb[0] as usize]];
        let mut dir = Direction::MATCH;

        if score < vgap_score {
            score = vgap_score;
            dir = Direction::gap_a(0 - vgap_pos);
        }

        if score < hgap_score[0] {
            score = hgap_score[0];
            dir = Direction::gap_b(col as i32 - hgap_pos[0]);
        }

        if score < 0 {
            curr_score[0] = 0;
            score = 0;
            dir_matrix[col * sb_len] = Direction::STOP;
        } else {
            curr_score[0] = score;
            dir_matrix[col * sb_len] = dir;
            if score >= max_score {
                max_row = 0;
                max_col = col;
                max_score = score;
            }
        }

        if score + gap_open >= vgap_score + gap_extend {
            vgap_score = score + gap_open;
            vgap_pos = 0;
        } else {
            vgap_score += gap_extend;
        }

        if score + gap_open >= hgap_score[0] + gap_extend {
            hgap_score[0] = score + gap_open;
            hgap_pos[0] = col as i32;
        } else {
            hgap_score[0] += gap_extend;
        }

        for row in 1..sb_len {
            score = prev_score[row - 1] + score_matrix[[sa[col] as usize, sb[row] as usize]];
            dir = Direction::MATCH;

            if score < vgap_score {
                score = vgap_score;
                dir = Direction::gap_a((row as i32) - vgap_pos);
            }

            if score < hgap_score[row] {
                score = hgap_score[row];
                dir = Direction::gap_b(col as i32 - hgap_pos[row]);
            }

            if score < 0 {
                curr_score[row] = 0;
                score = 0;
                dir_matrix[col * sb_len + row] = Direction::STOP;
            } else {
                curr_score[row] = score;
                dir_matrix[col * sb_len + row] = dir;
                if score >= max_score {
                    max_row = row;
                    max_col = col;
                    max_score = score;
                }
            }

            if score + gap_open >= vgap_score + gap_extend {
                vgap_score = score + gap_open;
                vgap_pos = row as i32;
            } else {
                vgap_score += gap_extend;
            }

            if score + gap_open >= hgap_score[row] + gap_extend {
                hgap_score[row] = score + gap_open;
                hgap_pos[row] = col as i32;
            } else {
                hgap_score[row] += gap_extend;
            }
        }
        std::mem::swap(&mut curr_score, &mut prev_score);
    }

    let align_frag = traceback(&dir_matrix, sb_len, max_col, max_row);

    let frag_count = align_frag.len();
    Ok(Alignment {
        align_frag: align_frag,
        frag_count: frag_count,
        score: max_score,
    })
}

#[pyfunction]
fn global_align(
    seqa: &[u8],
    seqb: &[u8],
    score_matrix: PyReadonlyArray2<i32>,
    gap_open: i32,
    gap_extend: i32,
) -> PyResult<Alignment> {
    let sa = seqa;
    let sb = seqb;
    let score_matrix = score_matrix.as_array();

    let sa_len = sa.len();
    let sb_len = sb.len();

    if sa_len == 0 || sb_len == 0 {
        return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
            "Input sequences cannot be empty.",
        ));
    }

    let mut curr_score = vec![0; sb_len];
    let mut prev_score = vec![0; sb_len];
    let mut dir_matrix = vec![Direction::MATCH; sa_len * sb_len];
    let mut hgap_pos = vec![0; sb_len];
    let mut hgap_score = vec![0; sb_len];

    for row in 0..sb_len {
        prev_score[row] = gap_open + gap_extend * row as i32;
        hgap_pos[row] = -1;
        hgap_score[row] = prev_score[row] + gap_open;
    }

    for col in 0..sa_len {
        let mut vgap_pos = -1;
        let mut vgap_score = gap_open * 2 + gap_extend * col as i32;

        let mut score = score_matrix[[sa[col] as usize, sb[0] as usize]]
            + (if col > 0 {
                gap_open + gap_extend * (col as i32 - 1)
            } else {
                0
            });
        let mut dir = Direction::MATCH;

        if score < vgap_score {
            score = vgap_score;
            dir = Direction::gap_a(0 - vgap_pos);
        }

        if score < hgap_score[0] {
            score = hgap_score[0];
            dir = Direction::gap_b(col as i32 - hgap_pos[0]);
        }

        curr_score[0] = score;
        dir_matrix[col * sb_len] = dir;

        if score + gap_open >= vgap_score + gap_extend {
            vgap_score = score + gap_open;
            vgap_pos = 0;
        } else {
            vgap_score += gap_extend;
        }

        if score + gap_open >= hgap_score[0] + gap_extend {
            hgap_score[0] = score + gap_open;
            hgap_pos[0] = col as i32;
        } else {
            hgap_score[0] += gap_extend;
        }

        for row in 1..sb_len {
            score = prev_score[row - 1] + score_matrix[[sa[col] as usize, sb[row] as usize]];
            dir = Direction::MATCH;

            if score < vgap_score {
                score = vgap_score;
                dir = Direction::gap_a((row as i32) - vgap_pos);
            }

            if score < hgap_score[row] {
                score = hgap_score[row];
                dir = Direction::gap_b(col as i32 - hgap_pos[row]);
            }

            curr_score[row] = score;
            dir_matrix[col * sb_len + row] = dir;

            if score + gap_open >= vgap_score + gap_extend {
                vgap_score = score + gap_open;
                vgap_pos = row as i32;
            } else {
                vgap_score += gap_extend;
            }

            if score + gap_open >= hgap_score[row] + gap_extend {
                hgap_score[row] = score + gap_open;
                hgap_pos[row] = col as i32;
            } else {
                hgap_score[row] += gap_extend;
            }
        }
        std::mem::swap(&mut curr_score, &mut prev_score);
    }

    let final_score = prev_score[sb_len - 1];
    let align_frag = global_traceback(&dir_matrix, sa_len, sa_len - 1, sb_len - 1);

    let frag_count = align_frag.len();
    Ok(Alignment {
        align_frag: align_frag,
        frag_count: frag_count,
        score: final_score,
    })
}

#[pyfunction]
fn glocal_align(
    seqa: &[u8],
    seqb: &[u8],
    score_matrix: PyReadonlyArray2<i32>,
    gap_open: i32,
    gap_extend: i32,
) -> PyResult<Alignment> {
    let sa = seqa;
    let sb = seqb;
    let score_matrix = score_matrix.as_array();

    let sa_len = sa.len();
    let sb_len = sb.len();

    if sa_len == 0 || sb_len == 0 {
        return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
            "Input sequences cannot be empty.",
        ));
    }

    let mut curr_score = vec![0; sb_len];
    let mut prev_score = vec![0; sb_len];
    let mut dir_matrix = vec![Direction::MATCH; sa_len * sb_len];
    let mut hgap_pos = vec![0; sb_len];
    let mut hgap_score = vec![0; sb_len];

    let mut max_score = 0;
    let mut max_row = 0;
    let mut max_col = 0;

    for row in 0..sb_len {
        prev_score[row] = std::i32::MIN / 2;
        hgap_pos[row] = -1;
        hgap_score[row] = std::i32::MIN / 2;
    }

    for col in 0..sa_len {
        let mut score = score_matrix[[sa[col] as usize, sb[0] as usize]];
        let dir = Direction::MATCH;

        curr_score[0] = score;
        dir_matrix[col * sb_len] = dir;

        let mut vgap_score = score + gap_open;
        let mut vgap_pos = 0;

        for row in 1..sb_len - 1 {
            score = prev_score[row - 1] + score_matrix[[sa[col] as usize, sb[row] as usize]];
            let mut dir = Direction::MATCH;

            if score < vgap_score {
                score = vgap_score;
                dir = Direction::gap_a((row as i32) - vgap_pos);
            }

            if score < hgap_score[row] {
                score = hgap_score[row];
                dir = Direction::gap_b(col as i32 - hgap_pos[row]);
            }

            curr_score[row] = score;
            dir_matrix[col * sb_len + row] = dir;

            if score + gap_open >= vgap_score + gap_extend {
                vgap_score = score + gap_open;
                vgap_pos = row as i32;
            } else {
                vgap_score += gap_extend;
            }

            if score + gap_open >= hgap_score[row] + gap_extend {
                hgap_score[row] = score + gap_open;
                hgap_pos[row] = col as i32;
            } else {
                hgap_score[row] += gap_extend;
            }
        }

        score = prev_score[sb_len - 2] + score_matrix[[sa[col] as usize, sb[sb_len - 1] as usize]];
        let dir = Direction::MATCH;

        curr_score[sb_len - 1] = score;
        dir_matrix[col * sb_len + sb_len - 1] = dir;

        if score >= max_score {
            max_row = sb_len - 1;
            max_col = col;
            max_score = score;
        }
        std::mem::swap(&mut curr_score, &mut prev_score);
    }
    let align_frag = traceback(&dir_matrix, sb_len, max_col, max_row);

    let frag_count = align_frag.len();
    Ok(Alignment {
        align_frag: align_frag,
        frag_count: frag_count,
        score: max_score,
    })
}

#[pyfunction]
fn overlap_align(
    seqa: &[u8],
    seqb: &[u8],
    score_matrix: PyReadonlyArray2<i32>,
    gap_open: i32,
    gap_extend: i32,
) -> PyResult<Alignment> {
    let sa = seqa;
    let sb = seqb;
    let score_matrix = score_matrix.as_array();

    let sa_len = sa.len();
    let sb_len = sb.len();

    if sa_len == 0 || sb_len == 0 {
        return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
            "Input sequences cannot be empty.",
        ));
    }

    let mut curr_score = vec![0; sb_len];
    let mut prev_score = vec![0; sb_len];
    let mut dir_matrix = vec![Direction::MATCH; sa_len * sb_len];
    let mut hgap_pos = vec![0; sb_len];
    let mut hgap_score = vec![0; sb_len];

    let mut max_score = std::i32::MIN;
    let mut max_row = -1;
    let mut max_col = -1;

    for row in 0..sb_len {
        curr_score[row] = score_matrix[[sa[0] as usize, sb[row] as usize]];
        hgap_pos[row] = 0;
        hgap_score[row] = curr_score[row] + gap_open;
        dir_matrix[row] = Direction::MATCH;
    }

    if curr_score[sb_len - 1] >= max_score {
        max_row = sb_len as i32 - 1;
        max_col = 0;
        max_score = curr_score[sb_len - 1];
    }

    std::mem::swap(&mut curr_score, &mut prev_score);

    for col in 1..sa_len - 1 {
        let mut score;
        let score_val = score_matrix[[sa[col] as usize, sb[0] as usize]];
        curr_score[0] = score_val;
        score = score_val;
        dir_matrix[col * sb_len] = Direction::MATCH;
        let mut vgap_pos = 0;
        let mut vgap_score = score + gap_open;

        for row in 1..sb_len {
            score = prev_score[row - 1] + score_matrix[[sa[col] as usize, sb[row] as usize]];
            let mut dir = Direction::MATCH;

            if score < vgap_score {
                score = vgap_score;
                dir = Direction::gap_a((row as i32) - vgap_pos);
            }

            if score < hgap_score[row] {
                score = hgap_score[row];
                dir = Direction::gap_b(col as i32 - hgap_pos[row]);
            }

            curr_score[row] = score;
            dir_matrix[col * sb_len + row] = dir;

            let is_gap_a = matches!(dir.kind(), DirectionKind::GapA(_));

            if !is_gap_a && score + gap_open >= vgap_score + gap_extend {
                vgap_score = score + gap_open;
                vgap_pos = row as i32;
            }
            else {
                vgap_score += gap_extend;
            }

            let is_gap_b =  matches!(dir.kind(), DirectionKind::GapB(_));

            if !is_gap_b && score + gap_open >= hgap_score[row] + gap_extend {
                hgap_score[row] = score + gap_open;
                hgap_pos[row] = col as i32;
            }
            else {
                hgap_score[row] += gap_extend;
            }
        }

        if curr_score[sb_len - 1] >= max_score {
            max_row = sb_len as i32 - 1;
            max_col = col as i32;
            max_score = curr_score[sb_len - 1];
        }
        std::mem::swap(&mut curr_score, &mut prev_score);
    }

    let col = sa_len - 1;
    let mut score = score_matrix[[sa[col] as usize, sb[0] as usize]];
    dir_matrix[col * sb_len] = Direction::MATCH;
    let mut vgap_pos = 0;
    let mut vgap_score = score + gap_open;

    if score >= max_score {
        max_row = 0;
        max_col = sa_len as i32 - 1;
        max_score = score;
    }

    for row in 1..sb_len {
        let mut dir = Direction::MATCH;
        score = prev_score[row - 1] + score_matrix[[sa[col] as usize, sb[row] as usize]];

        if score < vgap_score {
            score = vgap_score;
            dir = Direction::gap_a((row as i32) - vgap_pos);
        }

        if score < hgap_score[row] {
            score = hgap_score[row];
            dir = Direction::gap_b(col as i32 - hgap_pos[row]);
        }

        curr_score[row] = score;
        dir_matrix[col * sb_len + row] = dir;

        if score >= max_score {
            max_row = row as i32;
            max_col = sa_len as i32 - 1;
            max_score = score;
        }

        if score + gap_open >= vgap_score + gap_extend {
            vgap_score = score + gap_open;
            vgap_pos = row as i32;
        } else {
            vgap_score += gap_extend;
        }

        if score + gap_open >= hgap_score[row] + gap_extend {
            hgap_score[row] = score + gap_open;
            hgap_pos[row] = col as i32;
        } else {
            hgap_score[row] += gap_extend;
        }
    }

    if max_score == std::i32::MIN {
        return Ok(Alignment {
            align_frag: Vec::new(),
            frag_count: 0,
            score: 0,
        });
    }

    let final_max_col = if max_col >= 0 { max_col as usize } else { 0 };
    let final_max_row = if max_row >= 0 { max_row as usize } else { 0 };

    let align_frag = traceback(
        &dir_matrix,
        sb_len,
        final_max_col,
        final_max_row,
    );

    let frag_count = align_frag.len();
    Ok(Alignment {
        align_frag: align_frag,
        frag_count: frag_count,
        score: max_score,
    })
}

#[pymodule]
fn seq_align(_py: Python, m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(local_align, m)?)?;
    m.add_function(wrap_pyfunction!(global_align, m)?)?;
    m.add_function(wrap_pyfunction!(glocal_align, m)?)?;
    m.add_function(wrap_pyfunction!(overlap_align, m)?)?;
    m.add_class::<Alignment>()?;
    m.add_class::<AlignFrag>()?;
    m.add_class::<FragType>()?;
    Ok(())
}
