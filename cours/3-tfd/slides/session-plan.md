# TFD — déroulé de séance (2×1h)

## Why this document exists

Past chapters were built slide-first: dense decks that walk through the
poly section by section. Students struggle with that in a 1h slot —
there isn't time to both present and let ideas land. This chapter tries
something different: **plan the session, not the deck.** The deck is
whatever the story needs, no more. If a beat doesn't need a slide, it
doesn't get one (board, Ableton, Audacity, or just talking are all
valid "supports").

Two things are fixed points, everything else is negotiable:
- **One idea per beat**, stated as a single sentence a student could
  repeat back. Everything in that beat serves that one sentence.
- **The E-D-A loop** (Exemple → Définition → Application) drives beats
  2–4 identically, so students learn the *shape* of the session after
  beat 2 and can ride it through beat 4 without re-orienting each time.

**Sharpening rule (still in force):** every slide is exactly one
equation, one figure, or one definition — nothing shares a slide unless
the two things genuinely are one idea (flagged explicitly where that
happens). Beat 1 no longer respects the *slide count* this implied
(it's grown to 7 slides, see below) but every individual slide inside
it still respects the *density* rule — richer beat, still bare slides.

**This chapter now runs across two sessions**, split right after the
discret⟺périodique duality is established. The deck (`slides.tex`) has
been built to this plan and compiles clean — 29 pages: title, 6 section
dividers, 22 content frames.

## Session split

| Session | Covers | Ends on |
|---|---|---|
| **1** | Rappel → Fonction discrète → Discret ⟺ périodique → Bilan | The duality theorem, both directions, established — closed with its own conclusion |
| **2** | La TFD → Conclusion | The TFD itself, built as the duality applied twice |

## Timing overview — Session 1

Budget: ~50 min of content, slack folded into beat 3 (the one that
historically overruns via the board derivation).

| Beat | Target | Slides | If short on time |
|------|--------|--------|-------------------|
| 1. Rappel | 12 min | 7 | Cut the Dirac and convolution frames to one line each spoken instead of shown — they're recalled again in beat 3 anyway |
| 2. Fonction discrète | 8 min | 3 | Skip the Audacity capture slide, describe it instead |
| 3. Discret ⟺ périodique | 22 min | 4 | State the Dirac comb TF result with a one-line justification instead of the full board derivation |
| 4. Bilan (séance 1) | 6 min | 3 | Never cut — it's session 1's landing |

If genuinely behind after the board derivation, it's beat 1's newer
frames (Dirac, convolution) that get compressed to spoken asides, not
beat 4 — the bilan is the only thing in session 1 that must survive
intact, same logic as the final conclusion in session 2.

## Timing overview — Session 2

Budget: ~30 min of content — noticeably shorter than session 1. Flagged
as an open question below: this probably wants a recap opener that
doesn't exist yet.

| Beat | Target | Slides | If short on time |
|------|--------|--------|-------------------|
| (open: recap of session 1?) | — | 0 | See open questions |
| 5. La TFD | 20 min | 5 | Cut "why it matters" to a spoken line, drop the FFT aside entirely |
| 6. Conclusion | 5 min | 1 | Never cut — it's the payoff and the memory hook |

---

## The running example: one sound, three appearances

**One single audio artifact** is built once and reappears three times,
so students are looking at the same thing throughout instead of
resetting context each time:

- **Beat 1 (session 1)**: heard live, spectrum shown live — establishes
  "spectrum" viscerally, before any math.
- **Beat 2 (session 1)**: seen zoomed to individual samples — the same
  clip is now visual evidence for "discrete," not a new example.
- **Beat 6 (session 2)**: same clip's TFD, clean vs. leaked — the
  chapter's actual payoff made audible/visible on the same artifact one
  last time, a full session later.

The session gap (beats 1–2 in session 1, payoff in session 2) makes
re-using the *exact same clip* even more valuable than originally
planned — it's a concrete anchor students can be reminded of across the
week's gap ("remember that sound from last week?").

**Tool split — Ableton for live, Audacity for precise:**
- **Ableton Live**, beat 1 only: a track with a synth (Operator,
  sine-wave oscillator) + the built-in **Spectrum** device
  (Audio Effects → Spectrum) in analyzer view. Play a single held note
  live — one sharp peak appears in real time. Then play a chord —
  several peaks appear at once. This is a *performance*, not a
  rendered video: set it up once before class, then just play the
  keyboard live. No slide needed — Live's own window is the visual.
- **Audacity**, beats 2 and 6: same clip, rendered once (export the
  Ableton take, or generate directly in Audacity with Generate → Tone).
  - Beat 2: zoom in on the waveform until individual samples show as
    dots.
  - Beat 6: **Analyze → Plot Spectrum**, rectangular window, two
    prepared selections of the same tone — one an integer number of
    periods (clean peak), one not (leaked peak). Compute the exact
    selection lengths ahead of time, same verify-before-delivering
    discipline as the poly figures and TP3's parameters — still an
    open question below.

---

## SESSION 1

### Beat 1 — Rappel (now richer, 7 slides)

**One sentence:** *Everything so far — series, transform, Dirac,
convolution — has been building toward one move: decomposing signals
into sines and cosines, and today's chapter reuses every piece of it.*

Expanded from a single graph slide into an actual toolbox recap,
because beat 3's derivation reuses all of it and shouldn't have to
introduce it cold. Each slide is still bare — one line, one formula, or
one figure — the density rule didn't change, only the beat's scope.

1. **"On décompose"** — one sentence, no equation: *décomposer un
   signal comme une somme de sinus et de cosinus.* Sets the whole
   chapter's arc (series → TF → discret/périodique → TFD) as instances
   of one idea, before any of them are named.
2. **Graph** — Séries de Fourier → Transformée de Fourier → Discret ⟺
   Périodique → TFD. Same figure as before.
3. **Séries de Fourier** — tag "fonctions périodiques (période $T$)",
   then the complex-series formula and its coefficient formula
   together (one slide, tightly coupled — a series and its own
   coefficient formula is one concept, not two).
4. **Transformée de Fourier** — tag "fonctions non périodiques", then
   $X(f)=\int x(t)e^{-2j\pi ft}\,\mathrm dt$ alone. No inverse shown —
   not load-bearing today, mention it only if asked.
5. **L'impulsion de Dirac** — one figure (a smooth curve, a Dirac arrow
   at $t_0$, a dashed line to the picked-out value) plus the sifting
   formula, $\int x(t)\delta(t-t_0)\,\mathrm dt = x(t_0)$. This is the
   visual that makes "Dirac extracts one value" concrete before it gets
   reused in beat 3's board derivation.
6. **Produit de convolution** — the definition alone,
   $(x*y)(t)=\int x(\tau)y(t-\tau)\,\mathrm d\tau$.
7. **Dualité produit / convolution** — both directions,
   $x(t)y(t)\tf(X*Y)(f)$ and $(x*y)(t)\tf X(f)Y(f)$ — this is the exact
   tool beat 3's theorem proof needs, now sitting one beat away instead
   of being recalled from memory mid-derivation.

**Live, no slide:** the Ableton hook — under a minute, no commentary
beyond "this is what a spectrum *is*."

**The reframe, spoken, not on a slide:** Séries de Fourier was already
a special case of today's theorem, part b) — *periodic ⟹ discrete
spectrum* — proved without the general framing.

**Transition into beat 2:** "That sound was continuous. But the only
way a computer ever touches it is discretized. What *is* a discrete
signal, formally?"

### Beat 2 — Fonction discrète (unchanged, 3 slides)

**One sentence:** *A discrete signal is a function on $\Z$, full stop —
sampling is one origin among many, not the definition.*

1. **Exemple** — the Audacity zoomed-sample capture.
2. **Définition** — $x:\Z\to\R,\ n\mapsto x[n]$, alone.
3. **Application** — $u[n]$ and $\delta[n]$ side by side (the pairing
   exception: two canonical building blocks, contrasted).

**Transition into beat 3:** "So $x[n]$ is just a sequence. But we
already have a machine — the Fourier transform — built for continuous
signals. How do we connect the two?"

### Beat 3 — Discret ⟺ périodique (now 4 slides, temporal plots added)

**One sentence:** *Multiplying by a Dirac comb is how "discrete" enters
the continuous-Fourier world — and that one move is reversible.*

1. **Peigne de Dirac** — the comb figure (arrows), plus the sum formula
   $\Sha_{T_e}(t)=\sum_n\delta(t-nT_e)$.
2. **Board, live, no slide:** $\mathcal F[\Sha_{T_e}]$, full derivation
   — reusing beat 1's Dirac sifting property and the exponential↔Dirac
   pair, recalled out loud at the point each is used.
3. **Un résultat à connaître** — the boxed result,
   $\Sha_{T_e}(t)\tf F_e\Sha_{F_e}(f)$, revealed right after the board
   work concludes.
4. **Two frames, not one — this is the sharpening from your last
   round:** each direction of the theorem now gets its own frame with
   a temporal *and* a frequency panel side by side, not frequency
   alone:
   - **"$x$ discret ⟹ $X_e$ périodique"**: left panel, a dashed
     continuous envelope with solid sampling stems (visually: $x(t)
     \times \Sha_{T_e}(t)$) — right panel, the periodized-bump
     spectrum. The time panel makes "discretizing" concrete instead of
     leaving it as a frequency-only abstraction.
   - **"$x$ périodique ⟹ $X$ discret"**: left panel, an actual periodic
     wave — right panel, the Dirac-comb spectrum.

   Splitting the old single two-panel frame into two frames (one per
   direction) was necessary once each direction got its own time+freq
   pair — cramming four panels onto one slide would have broken the
   density rule outright.

**Gut-check, spoken, no slide:** "so if I *sample* the sound, what
happens to its spectrum?" (periodizes it).

**Transition into beat 4:** "We now know discrete ⟺ periodic in
general — in both directions, in time and in frequency."

### Beat 4 — Bilan (séance 1) — NEW

**One sentence:** *Discret ⟺ périodique is proved, both directions —
next time we apply it a second time and get the TFD for free.*

This beat didn't exist in the previous draft; it's the direct answer
to splitting the chapter into two sessions. Its only job is to close
session 1 as a complete unit, not a cliffhanger mid-derivation.

1. **"Où on en est"** — the graph again, but now with the *third* node
   (Discret ⟺ Périodique) filled instead of TFD — visually: "we got
   this far today, not further." One line underneath: *discret ⟺
   périodique, établi dans les deux sens.*
2. **"Même théorème, autre paire"** — added per your remark that the
   $(n,\nu)$ change of variable needed to land clearly before the
   week's gap, not just get mentioned in passing. One correspondence,
   labeled, nothing else: $t\leftrightarrow n$, $f\leftrightarrow\nu$,
   each pair captioned in three words (temps continu/indice discret,
   fréquence/fréquence normalisée), then one closing line — *le
   théorème de dualité s'applique tel quel.* This is the actual bridge
   to session 2: without it, "on applique la dualité une seconde fois"
   is just an assertion; with it, students have already seen which
   variables get swapped.
3. **"La semaine prochaine"** — now just the one-line teaser, *on
   applique cette dualité une seconde fois pour construire la TFD* —
   the $(n,\nu)$ detail moved to slide 2 above, so this one stays pure
   teaser without repeating it.

**Do not cut any of these three slides** — same rule as the final
conclusion: this beat's entire job is making session 1 feel finished,
not truncated, and slide 2 specifically is load-bearing for session 2's
opening move.

---

## SESSION 2

### Beat 5 — La Transformée de Fourier Discrète (unchanged, 5 slides)

**One sentence:** *A finite signal only gets a finite spectrum if we're
willing to imagine it repeating forever — that's the one trick behind
the TFD.*

1. **"?"** — full-bleed question, no equation: *le théorème exige
   périodique. on a fini. et maintenant ?* Let it sit before answering
   out loud: we *choose* to treat the finite sequence as one period of
   an implicit periodic signal, purely so beat 3's theorem applies.
2. **Fréquence normalisée** — both lines together (one concept: a
   definition and its immediate periodicity consequence),
   $\nu=f/F_e$ and $X(\nu)=\sum x[n]e^{-2j\pi\nu n}$ (période 1).
3. **TFD** — $X[k]=\sum_{n=0}^{N-1}x[n]e^{-2j\pi kn/N}$ alone.
4. **TFD inverse** — its own slide, sequential reveal after the direct
   form, not shared with it.
5. **Pourquoi c'est important** — one figure (spectrum-analyzer /
   oscilloscope FFT-mode capture), one line: *c'est ce calcul,
   partout.* FFT complexity remark spoken only, first cut if short.

**Transition into beat 6:** "So — back to that sound clip."

### Beat 6 — Conclusion (unchanged, 1 slide)

**One sentence:** *Discrete ⟺ periodic, applied twice, is the entire
chapter.*

**Live, no slide:** the Audacity payoff — same tone, clean vs. leaked
spectrum, in Audacity's own window.

**Slide (shown last, left on screen while class packs up):** the graph
with TFD now filled in, plus: *discret ⟺ périodique, appliqué deux
fois : la TFD.*

**Never cut** — it's the only beat whose entire job is making the
chapter memorable rather than correct.

---

## Open questions (need your call)

1. **Session 2 opener:** right now session 2 starts cold on "?" (beat
   5, slide 1) with no recap of session 1's duality theorem after a
   week's gap. Do you want a short recap slide (e.g. the beat-4 graph
   again, or a one-line restatement of the theorem) at the top of
   session 2, mirroring beat 1's own "previously on ch1/ch2" framing? I
   didn't add one unprompted since you gave a specific beat list, but
   session 2's ~30 min budget has room for it.
2. **Board derivation time for the Dirac comb TF (beat 3):** how long
   does this actually take you live? The 22 min budget assumes a good
   chunk goes to it — if it runs longer, it eats into beat 1's newer
   frames (Dirac/convolution get compressed to spoken asides first).
3. **Audacity dry run:** do you want me to work out the actual tone
   frequency / sample rate / selection-length numbers for the beat-6
   clean-vs-leaked demo now (same numeric-verification pass as TP3's
   parameters), or set that up yourself once the deck exists?

## Explicitly out of scope for this session

- Shannon/Nyquist, ideal reconstruction, aliasing by name — already
  covered thoroughly in elec-s3; this chapter stays on the pure
  discret⟺périodique duality and doesn't re-teach sampling theory.
- Windowing as a separate topic — fuite spectrale is explained here via
  the implicit-periodic-extension route (beat 5/6), not via
  window-function convolution (elec-s3's route) — don't let the two
  explanations blur together in the oral delivery.
