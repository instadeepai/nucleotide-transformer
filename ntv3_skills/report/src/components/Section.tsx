import type { ReactNode } from "react";

type SectionProps = {
  id: string;
  title: string;
  children: ReactNode;
};

export function Section({ id, title, children }: SectionProps) {
  return (
    <section className="section" id={id}>
      <header className="section__head">
        <h2>{title}</h2>
      </header>
      <div className="section__body">{children}</div>
    </section>
  );
}
