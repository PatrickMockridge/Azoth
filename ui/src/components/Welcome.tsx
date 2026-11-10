/**
 * The first visit, over the demo rather than instead of it.
 *
 * **A card and not a modal.** The demo is already behind it, so the page a newcomer lands on is
 * a working flowsheet with a short explanation, not a gate that hides the thing they came to
 * see. Every action dismisses the card, and the browser remembers that it was dealt with.
 */

import type { Example } from "../examples";

export interface WelcomeProps {
  examples: Example[];
  onDemo: () => void;
  onBlank: () => void;
  onExample: (document: string) => void;
  onHelp: () => void;
  onDismiss: () => void;
}

export function Welcome({ examples, onDemo, onBlank, onExample, onHelp, onDismiss }: WelcomeProps) {
  return (
    <section className="welcome" aria-labelledby="welcome-title">
      <button
        type="button"
        className="welcome-close"
        aria-label="dismiss the welcome card"
        onClick={onDismiss}
      >
        ×
      </button>
      <h2 id="welcome-title">A flowsheet in the browser</h2>
      <p>
        The demo behind this card is a feed, a pump, a heater, a separator and a liquid recycle.
        Press <strong>Solve</strong> to run it, then click a unit operation or a stream to read
        its values.
      </p>
      <div className="welcome-actions">
        <button type="button" className="primary" onClick={onDemo}>
          Open the demo
        </button>
        <button type="button" onClick={onBlank}>
          Start blank
        </button>
        <button type="button" onClick={onHelp}>
          Read the 2-minute guide
        </button>
      </div>
      {examples.length === 0 ? null : (
        <>
          <h3>Or start from a worked case</h3>
          <div className="welcome-examples">
            {examples.map((example) => (
              <button
                key={example.file}
                type="button"
                onClick={() => onExample(example.document)}
                title={example.description}
              >
                {example.title}
              </button>
            ))}
          </div>
        </>
      )}
      <p className="note">
        The port is incomplete and every model is marked unverified or partially verified. These
        are screening values, not design deliverables.
      </p>
    </section>
  );
}
