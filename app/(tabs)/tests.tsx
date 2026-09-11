import { FlatList, StyleSheet, Text, View } from 'react-native';
import { ScreenContainer } from '@/components/screen-container';
import { IconSymbol } from '@/components/ui/icon-symbol';

type Contract = { id: string; version: string; engine: string; evidence: string; state: string; checks: string };
const contracts: Contract[] = [
  { id: 'R4.1-QUALIFICATION', version: 'v1.0.0', engine: 'OpenCascade adapter', evidence: 'ARTIFACT_PROOF', state: 'REGISTERED', checks: 'Identity · SHA · STEP · B-Rep · Topology · Gate' },
  { id: 'STEP-TOPOLOGY', version: 'v0.1.0', engine: 'OpenCascade adapter', evidence: 'GEOMETRY_PROOF', state: 'REGISTERED', checks: 'Schema · Bodies · Solids · Topology · Coverage' },
  { id: 'STRUCTURAL-FUTURE', version: 'v0.0.0', engine: 'OpenRadioss adapter', evidence: 'CAE_PROOF', state: 'NOT_PROVEN', checks: 'No solver runtime or engineering input available' },
];

export default function TestsScreen() {
  return <ScreenContainer containerClassName="bg-[#08111F]" className="px-5" edges={['top', 'left', 'right']}>
    <FlatList data={contracts} keyExtractor={(item) => item.id} contentContainerStyle={styles.content} ListHeaderComponent={<View style={styles.header}><Text style={styles.eyebrow}>CONTRACT REGISTRY</Text><Text style={styles.title}>Tests</Text><Text style={styles.subtitle}>A contract is not a result. Version, engine, checks, acceptance, evidence class, and gate policy stay explicit.</Text></View>} renderItem={({ item }) => <View style={styles.card}><View style={styles.cardTop}><View style={styles.icon}><IconSymbol name="checkmark.seal.fill" size={17} color="#2DE0B2" /></View><View style={styles.heading}><Text style={styles.id}>{item.id}</Text><Text style={styles.version}>{item.version} · {item.state}</Text></View></View><Text style={styles.engine}>ENGINE · {item.engine}</Text><Text style={styles.checks}>{item.checks}</Text><View style={styles.footer}><Text style={styles.evidence}>EVIDENCE · {item.evidence}</Text><Text style={styles.gate}>GATE POLICY · BLOCK UPSTREAM</Text></View></View>} />
  </ScreenContainer>;
}

const styles = StyleSheet.create({ content: { paddingTop: 16, paddingBottom: 38, gap: 12 }, header: { gap: 7, marginBottom: 4 }, eyebrow: { color: '#8DA0B8', fontSize: 9, fontWeight: '900', letterSpacing: 1.4 }, title: { color: '#F5F8FC', fontSize: 27, fontWeight: '900' }, subtitle: { color: '#8192A9', fontSize: 11, lineHeight: 17 }, card: { backgroundColor: '#101D31', borderColor: '#22334A', borderWidth: 1, borderRadius: 16, padding: 14, gap: 10 }, cardTop: { flexDirection: 'row', gap: 10, alignItems: 'center' }, icon: { width: 34, height: 34, borderRadius: 11, backgroundColor: '#173833', alignItems: 'center', justifyContent: 'center' }, heading: { flex: 1 }, id: { color: '#F5F8FC', fontSize: 13, fontWeight: '900' }, version: { color: '#2DE0B2', fontSize: 9, fontWeight: '900', marginTop: 3, letterSpacing: 0.6 }, engine: { color: '#FFBF69', fontSize: 9, fontWeight: '900', letterSpacing: 0.7 }, checks: { color: '#A0B0C2', fontSize: 11, lineHeight: 17 }, footer: { borderTopColor: '#22334A', borderTopWidth: 1, paddingTop: 9, gap: 5 }, evidence: { color: '#8DA0B8', fontSize: 8, fontWeight: '800', letterSpacing: 0.6 }, gate: { color: '#FF6B76', fontSize: 8, fontWeight: '900', letterSpacing: 0.6 },
});
