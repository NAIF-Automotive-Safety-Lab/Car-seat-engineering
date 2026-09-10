import { SymbolView, type SymbolViewProps } from 'expo-symbols';
import { type StyleProp, type ViewStyle } from 'react-native';

export type IconSymbolName = keyof typeof MAPPING;
const MAPPING = {
  'house.fill': 'house.fill', 'square.grid.2x2.fill': 'square.grid.2x2.fill', 'exclamationmark.triangle.fill': 'exclamationmark.triangle.fill', 'cpu.fill': 'cpu.fill', 'list.bullet.rectangle.portrait.fill': 'list.bullet.rectangle.portrait.fill', 'shield.lefthalf.filled': 'shield.lefthalf.filled', 'arrow.clockwise': 'arrow.clockwise', 'account-tree.fill': 'point.3.connected.trianglepath.dotted', 'arrow.triangle.2.circlepath': 'arrow.triangle.2.circlepath', 'waveform.path.ecg': 'waveform.path.ecg', 'checkmark.seal.fill': 'checkmark.seal.fill', 'doc.text.magnifyingglass': 'doc.text.magnifyingglass', 'chevron.right': 'chevron.right', 'lock.fill': 'lock.fill', 'hammer.fill': 'hammer.fill', 'checkmark': 'checkmark', 'clock.fill': 'clock.fill', 'person.crop.circle.fill': 'person.crop.circle.fill',
} as const;

export function IconSymbol({ name, size = 24, color, style, weight = 'regular' }: { name: IconSymbolName; size?: number; color: string; style?: StyleProp<ViewStyle>; weight?: SymbolViewProps['weight'] }) {
  return <SymbolView name={MAPPING[name] as SymbolViewProps['name']} weight={weight} tintColor={color} style={[{ width: size, height: size }, style]} />;
}
