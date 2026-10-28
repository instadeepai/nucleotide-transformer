import catalogFile from "./ucscCatalog.json" with { type: "json" };
import type { CatalogTrack } from "./ucscGroups";

type Catalog = {
  genome: string;
  nTop: number;
  tracks: CatalogTrack[];
};

export const UCSC_CATALOG = catalogFile as Catalog;

export function findCatalogTrack(trackId: string, genome: string): {
  id: string;
  label: string;
  type: string;
  group: string;
} | null {
  if (genome !== UCSC_CATALOG.genome) return null;
  for (const track of UCSC_CATALOG.tracks) {
    if (track.id === trackId) {
      return {
        id: track.id,
        label: track.label,
        type: track.type,
        group: track.group,
      };
    }
    const leaf = track.leaves.find((candidate) => candidate.id === trackId);
    if (leaf) {
      return {
        id: leaf.id,
        label: leaf.label,
        type: leaf.type,
        group: track.group,
      };
    }
  }
  return null;
}
