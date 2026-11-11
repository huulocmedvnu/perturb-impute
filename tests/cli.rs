//! End-to-end test of the perturb-impute command line on the toy data in tests/data:
//! panel selection (standard and low-budget), reference fitting, imputation with INFO and
//! out-of-reference scores, and projection of whole-transcriptome profiles on reference programs.
use std::collections::HashMap;
use std::path::{Path, PathBuf};
use std::process::Command;

fn run(args: &[&str]) {
    let out = Command::new(env!("CARGO_BIN_EXE_perturb-impute")).args(args).output().expect("run perturb-impute");
    assert!(out.status.success(), "perturb-impute {:?} failed: {}", args, String::from_utf8_lossy(&out.stderr));
}
fn data(name: &str) -> String { Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data").join(name).to_str().unwrap().to_string() }
fn s(p: &PathBuf) -> &str { p.to_str().unwrap() }
/// TSV matrix: (row names, column names, values by row)
fn matrix(path: &str) -> (Vec<String>, Vec<String>, Vec<Vec<f64>>) {
    let txt = std::fs::read_to_string(path).unwrap_or_else(|e| panic!("{path}: {e}"));
    let mut lines = txt.lines();
    let cols: Vec<String> = lines.next().unwrap().split('\t').skip(1).map(String::from).collect();
    let (mut rows, mut v) = (vec![], vec![]);
    for l in lines { let mut it = l.split('\t'); rows.push(it.next().unwrap().to_string()); v.push(it.map(|x| x.parse::<f64>().unwrap()).collect()); }
    (rows, cols, v)
}
fn pearson(x: &[f64], y: &[f64]) -> f64 {
    let n = x.len() as f64;
    let (mx, my) = (x.iter().sum::<f64>() / n, y.iter().sum::<f64>() / n);
    let (mut sxy, mut sxx, mut syy) = (0.0, 0.0, 0.0);
    for (a, b) in x.iter().zip(y) { sxy += (a - mx) * (b - my); sxx += (a - mx).powi(2); syy += (b - my).powi(2); }
    sxy / (sxx * syy).sqrt()
}
fn median(mut v: Vec<f64>) -> f64 { v.sort_by(|a, b| a.partial_cmp(b).unwrap()); v[v.len() / 2] }
/// per-perturbation correlation with the truth on the genes named in `cols` (the out-of-reference row excluded)
fn row_r(rows: &[String], cols: &[String], v: &[Vec<f64>], truth: &(Vec<String>, Vec<String>, Vec<Vec<f64>>)) -> Vec<f64> {
    let ti: HashMap<&str, usize> = truth.1.iter().enumerate().map(|(i, g)| (g.as_str(), i)).collect();
    let tr: HashMap<&str, usize> = truth.0.iter().enumerate().map(|(i, p)| (p.as_str(), i)).collect();
    rows.iter().zip(v).filter(|(p, _)| p.as_str() != "OUT_OF_REFERENCE").map(|(p, x)| {
        let t: Vec<f64> = cols.iter().map(|g| truth.2[tr[p.as_str()]][ti[g.as_str()]]).collect();
        pearson(x, &t)
    }).collect()
}

#[test]
fn select_fit_impute_project() {
    let dir = std::env::temp_dir().join(format!("perturb_impute_cli_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let (std_panel, low_panel, model, model_low, imp, proj) = (dir.join("standard.tsv"), dir.join("low_budget.tsv"), dir.join("model.bin"), dir.join("model_low.bin"), dir.join("imputed"), dir.join("projected.tsv"));
    let (reference, umi) = (data("toy_reference.tsv"), data("toy_umi_fraction.tsv"));

    // panel selection: greedy MMSE, and the UMI-cost variant
    run(&["select", "--ref", &reference, "--panel-size", "40", "--umi-fraction", &umi, "--out", s(&std_panel)]);
    run(&["select", "--ref", &reference, "--panel-size", "40", "--mode", "low-budget", "--umi-fraction", &umi, "--out", s(&low_panel)]);
    let panel = |p: &PathBuf| -> (Vec<String>, f64) {
        let t = std::fs::read_to_string(p).unwrap();
        let rows: Vec<Vec<String>> = t.lines().skip(1).map(|l| l.split('\t').map(String::from).collect()).collect();
        (rows.iter().map(|r| r[1].clone()).collect(), rows.last().unwrap()[4].parse().unwrap())
    };
    let ((g_std, share_std), (g_low, share_low)) = (panel(&std_panel), panel(&low_panel));
    assert_eq!((g_std.len(), g_low.len()), (40, 40));
    assert!(share_low < share_std, "the UMI-cost panel should hold a smaller share of UMIs ({share_low} vs {share_std})");

    // reference models, with and without noise weights for a sequencing budget
    run(&["fit-reference", "--ref", &reference, "--panel", s(&std_panel), "--groups", &data("toy_groups.tsv"), "--out", s(&model)]);
    run(&["fit-reference", "--ref", &reference, "--panel", s(&low_panel), "--groups", &data("toy_groups.tsv"), "--noise", &data("toy_noise.tsv"), "--budget", "0.02", "--out", s(&model_low)]);
    assert!(std::fs::metadata(&model_low).unwrap().len() > 0);

    // imputation of the 260 unmeasured genes of 41 screen profiles
    run(&["impute", "--model", s(&model), "--input", &data("toy_screen_input.tsv"), "--out-dir", s(&imp)]);
    let truth = matrix(&data("toy_screen_truth.tsv"));
    let (rows, cols, z) = matrix(s(&imp.join("imputed_zscores.tsv")));
    assert_eq!((rows.len(), cols.len()), (41, 260));
    let r_imp = median(row_r(&rows, &cols, &z, &truth));
    assert!(r_imp > 0.9, "median correlation of imputed genes with the truth {r_imp}");
    let (_, _, lfsr) = matrix(s(&imp.join("imputed_lfsr.tsv")));
    assert!(lfsr.iter().flatten().all(|p| (0.0..=0.5).contains(p)), "local false sign rates lie in [0, 0.5]");
    let info: Vec<f64> = std::fs::read_to_string(imp.join("gene_info.tsv")).unwrap().lines().skip(1).map(|l| l.split('\t').nth(1).unwrap().parse().unwrap()).collect();
    assert!(median(info) > 0.9, "median INFO score");
    // the profile without program structure has the highest out-of-reference score
    let qc: Vec<(String, f64)> = std::fs::read_to_string(imp.join("perturbation_qc.tsv")).unwrap().lines().skip(1).map(|l| { let f: Vec<&str> = l.split('\t').collect(); (f[0].to_string(), f[1].parse().unwrap()) }).collect();
    let top = qc.iter().max_by(|a, b| a.1.partial_cmp(&b.1).unwrap()).unwrap();
    assert_eq!(top.0, "OUT_OF_REFERENCE");
    assert!(qc.iter().all(|q| (0.0..=1.0).contains(&q.1)));

    // projection on the 8 reference programs removes measurement noise
    run(&["project", "--ref", &reference, "--input", &data("toy_screen_input.tsv"), "--programs", "8", "--out", s(&proj)]);
    let (prow, pcol, pz) = matrix(s(&proj));
    let input = matrix(&data("toy_screen_input.tsv"));
    let (r_proj, r_raw) = (median(row_r(&prow, &pcol, &pz, &truth)), median(row_r(&input.0, &input.1, &input.2, &truth)));
    assert!(r_proj > r_raw, "projected profiles ({r_proj}) should agree with the truth better than measured ones ({r_raw})");
    std::fs::remove_dir_all(&dir).ok();
}

/// guide-target prior: reference knockdowns named after genes give the value of the target gene and of genes
/// with a TSS within 2 kb; `impute --guide-target` overwrites exactly these entries and nothing else
#[test]
fn guide_target_prior() {
    let dir = std::env::temp_dir().join(format!("perturb_impute_prior_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    // reference: knockdown i is named after gene G{i}, with an on-target effect of -6; gene i sits at position 5000 * i,
    // and the odd gene i + 1 at 1000 bp from the even gene i, so every even target has one neighbour within 2 kb
    let txt = std::fs::read_to_string(data("toy_reference.tsv")).unwrap();
    let mut lines = txt.lines();
    let header = lines.next().unwrap().to_string();
    let genes: Vec<&str> = header.split('\t').skip(1).collect();
    let mut out = vec![header.clone()];
    for (i, l) in lines.enumerate() {
        let mut f: Vec<String> = l.split('\t').map(String::from).collect();
        f[0] = genes[i].to_string();
        f[i + 1] = "-6.000".into();
        if i % 2 == 0 { f[i + 2] = "-2.500".into(); }
        out.push(f.join("\t"));
    }
    let (reference, tss, expr, panel, model, imp, imp0) = (dir.join("ref.tsv"), dir.join("tss.tsv"), dir.join("expr.tsv"), dir.join("panel.tsv"), dir.join("model.bin"), dir.join("with_prior"), dir.join("plain"));
    std::fs::write(&reference, out.join("\n") + "\n").unwrap();
    let pos = |i: usize| if i % 2 == 0 { 5000 * i } else { 5000 * (i - 1) + 1000 };
    std::fs::write(&tss, format!("gene\tchr\tstrand\ttss\n{}\n", genes.iter().enumerate().map(|(i, g)| format!("{g}\tchr1\t+\t{}", pos(i))).collect::<Vec<_>>().join("\n"))).unwrap();
    std::fs::write(&expr, genes.iter().enumerate().map(|(i, g)| format!("{g}\t{}", i + 1)).collect::<Vec<_>>().join("\n") + "\n").unwrap();
    run(&["select", "--ref", s(&reference), "--panel-size", "40", "--out", s(&panel)]);
    run(&["fit-reference", "--ref", s(&reference), "--panel", s(&panel), "--tss-annotation", s(&tss), "--expression", s(&expr), "--out", s(&model)]);
    // screen: knockdowns of genes that are not in the reference (indices 240 and above)
    let scr = std::fs::read_to_string(data("toy_screen_input.tsv")).unwrap();
    let rows: Vec<String> = scr.lines().enumerate().map(|(i, l)| if i == 0 { l.to_string() } else { let mut f: Vec<&str> = l.split('\t').collect(); f[0] = genes[240 + i - 1]; f.join("\t") }).collect();
    let screen = dir.join("screen.tsv");
    std::fs::write(&screen, rows.join("\n") + "\n").unwrap();
    run(&["impute", "--model", s(&model), "--input", s(&screen), "--out-dir", s(&imp0)]);
    run(&["impute", "--model", s(&model), "--input", s(&screen), "--guide-target", "--out-dir", s(&imp)]);
    let (r0, c0, z0) = matrix(s(&imp0.join("imputed_zscores.tsv")));
    let (_, _, z1) = matrix(s(&imp.join("imputed_zscores.tsv")));
    let panel_genes: Vec<String> = std::fs::read_to_string(&panel).unwrap().lines().skip(1).map(|l| l.split('\t').nth(1).unwrap().to_string()).collect();
    let (mut changed, mut expected) = (0, 0);
    for (i, name) in r0.iter().enumerate() {
        let t: usize = name[1..].parse().unwrap();
        for (j, gname) in c0.iter().enumerate() {
            let q: usize = gname[1..].parse().unwrap();
            let is_target = q == t;
            let is_near = q != t && (pos(q) as i64 - pos(t) as i64).abs() < 2000;
            if is_target || is_near {
                expected += 1;
                // on-target effects in the reference are all -6; neighbours of even targets are -2.5, of odd targets ordinary values
                if is_target { assert!((z1[i][j] + 6.0).abs() < 1e-3, "target {name}: {}", z1[i][j]); }
            } else { assert_eq!(z0[i][j], z1[i][j], "entry outside the prior changed"); }
            if z0[i][j] != z1[i][j] { changed += 1; }
        }
    }
    assert!(expected > 0 && changed > 0 && changed <= expected);
    let listed = std::fs::read_to_string(imp.join("prior_applied.tsv")).unwrap().lines().count() - 1;
    assert_eq!(listed, expected, "prior_applied.tsv lists every overwritten entry");
    assert!(!panel_genes.is_empty());
    // a model without the prior refuses the option
    let plain = dir.join("plain.bin");
    run(&["fit-reference", "--ref", s(&reference), "--panel", s(&panel), "--out", s(&plain)]);
    let o = Command::new(env!("CARGO_BIN_EXE_perturb-impute")).args(["impute", "--model", s(&plain), "--input", s(&screen), "--guide-target", "--out-dir", s(&imp)]).output().unwrap();
    assert!(!o.status.success());
    std::fs::remove_dir_all(&dir).ok();
}
