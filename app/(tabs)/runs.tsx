import { FlatList, StyleSheet, Text, View } from 'react-native';
import { ScreenContainer } from '@/components/screen-container';
import { IconSymbol } from '@/components/ui/icon-symbol';

type Run = { id: string; test: string; state: 'BLOCKED' | 'QUEUED' | 'NOT_PROVEN'; detail: string; gate: string };
const runs: Run[] = [
  { id: 'K-RUN-0001', test: 'R4.1-QUALIFICATION · v1.0.0', state: 'BLOCKED', detail: 'Canonical artifact SHA is unavailable.', gate: 'BLOCK BEFORE EXECUTION' },
  { id: 'K-RUN-0002', test: 'STEP-TOPOLOGY · v0.1.0', state: 'NOT_PROVEN', detail: 'Adapter is registered; no solver result exists.', gate: 'EVIDENCE REQUIRED' },
];

const tone: Record<Run['state'], string> = { BLOCKED: '#FF6B76', QUEUED: '#FFBF69', NOT_PROVEN: '#8DA0B8' };

export default function RunsScreen() {
  return <ScreenContainer containerClassName="bg-[#08111F]" className="px-5" edges={['top', 'left', 'right']}>
    <FlatList data={runs} keyExtractor={(item) => item.id} contentContainerStyle={styles.content} ListHeaderComponent={<View style={styles.header}><Text style={styles.eyebrow}>EXECUTION BUS</Text><Text style={styles.title}>Kernel runs</Text><Text style={styles.subtitle}>The UI can inspect runs; only the trusted kernel may execute them.</Text><View style={styles.policy}><IconSymbol name="shield.lefthalf.filled" size={18} color="#2DE0B2" /><View><Text style={styles.policyTitle}>No silent skip</Text><Text style={styles.policyCopy}>Every blocked or unverified state remains visible.</Text></View></View><View style={styles.legend}><Text style={styles.legendText}>REGISTERED</Text><Text style={styles.legendText}>EXECUTED</Text><Text style={styles.legendText}>VERIFIED</Text><Text style={styles.legendText}>PASS</Text></View></View>} renderItem={({ item }) => <View style={styles.card}><View style={styles.cardTop}><Text style={styles.id}>{item.id}</Text><Text style={[styles.state, { color: tone[item.state] }]}>{item.state}</Text></View><Text style={styles.test}>{item.test}</Text><Text style={styles.detail}>{item.detail}</Text><View style={styles.divider} /><Text style={styles.gate}>{item.gate}</Text></View>} />
  </ScreenContainer>;
}

const styles = StyleSheet.create({ content: { paddingTop: 16, paddingBottom: 38, gap: 12 }, header: { gap: 7, marginBottom: 4 }, eyebrow: { color: '#8DA0B8', fontSize: 9, fontWeight: '900', letterSpacing: 1.4 }, title: { color: '#F5F8FC', fontSize: 27, fontWeight: '900' }, subtitle: { color: '#8192A9', fontSize: 11, lineHeight: 17 }, policy: { flexDirection: 'row', gap: 11, alignItems: 'center', backgroundColor: '#12352F', borderColor: '#1E6C5B', borderWidth: 1, borderRadius: 15, padding: 13, marginTop: 8 }, policyTitle: { color: '#D9FFF5', fontWeight: '900', fontSize: 12 }, policyCopy: { color: '#8DC9BC', fontSize: 10, marginTop: 3 }, legend: { flexDirection: 'row', justifyContent: 'space-between', marginTop: 8 }, legendText: { color: '#667991', fontSize: 8, fontWeight: '800', letterSpacing: 0.6 }, card: { backgroundColor: '#101D31', borderColor: '#22334A', borderWidth: 1, borderRadius: 16, padding: 14, gap: 7 }, cardTop: { flexDirection: 'row', justifyContent: 'space-between' }, id: { color: '#8DA0B8', fontSize: 9, fontWeight: '900', letterSpacing: 1 }, state: { fontSize: 9, fontWeight: '900', letterSpacing: 0.7 }, test: { color: '#F5F8FC', fontSize: 14, fontWeight: '900' }, detail: { color: '#A0B0C2', fontSize: 11, lineHeight: 17 }, divider: { height: 1, backgroundColor: '#22334A', marginTop: 3 }, gate: { color: '#FFBF69', fontSize: 9, fontWeight: '900', letterSpacing: 0.7 },
});
