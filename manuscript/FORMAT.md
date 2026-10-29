# Manuscript source format (converted to the narrated reader AND compiled to PDF)

Language: English only. Style: Codex skill ~/.codex/skills/mathematical-textbook/SKILL.md (read it) — numbered Definition / Proposition / Theorem / Lemma / Corollary / Example / Remark; Proof and Example in their own blocks; formula chains carry the argument; short connecting prose. The user explicitly asked that the STORY be clear: each chapter opens with 2-4 sentences saying what problem it solves and why it is not trivial, and the overview chapter tells the whole story. No defensive wording ("only", "does not", "not a ..."); state objects and results affirmatively. Source notation wins over the skill's abstract-index default (user's explicit request). Every symbol: from sources/ledger_SYNTHESIS.md (cite the source the first time, e.g. "[C07, Lemma 10.2.1]") or written out in full. No audit notes, script names, residual counts or file paths in the manuscript.

File: one chapter per file manuscript/chNN_slug.tex. Allowed markup ONLY:
- First line: \chapter{Title}
- \section{Title}, \subsection{Title}
- Paragraphs separated by a blank line; inline math $...$ (no $$).
- Display math: \begin{equation} ... \label{eq:...} \end{equation} or \[ ... \]. One display per environment; multi-line via \begin{aligned} ... \end{aligned} inside. KaTeX-compatible LaTeX only (amsmath core, \operatorname, \mathrm, \mathfrak, \mathscr, \mathcal, \not, \partial, \hbar, \varepsilon, \langle\rangle, \underbrace, \overbrace, cases, pmatrix). No custom macros, no \slashed, no \text with math inside braces of other macros.
- Theorem-like: \begin{definition}[Title] ... \end{definition}; same for proposition, theorem, lemma, corollary, example, remark, axiom (= external input, give exact citation), proof. Labels: \label{def:...}, \label{prop:...}; references \ref{...}, \eqref{...}.
- Figures: \fig{diagrams/composite/<file>.png}{Caption text with $math$.} on its own line.
- Tables: \begin{tabular}{l|c|c} ... \end{tabular} with &, \\, \hline; cells may contain $math$; optional caption line right after: \tabcaption{...}.
- Lists: \begin{itemize} \item ... \end{itemize} (no nesting).
- Citations: plain text in square brackets, e.g. [ESW, (31)].
Nothing else (no \cite, no \footnote, no \newcommand, no \input).
