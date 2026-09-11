import MaterialIcons from '@expo/vector-icons/MaterialIcons';
import { ComponentProps } from 'react';
import { OpaqueColorValue, type StyleProp, type TextStyle } from 'react-native';

type IconName = ComponentProps<typeof MaterialIcons>['name'];
export type IconSymbolName = keyof typeof MAPPING;

const MAPPING = {
  'house.fill': 'home', 'square.grid.2x2.fill': 'dashboard', 'exclamationmark.triangle.fill': 'warning', 'cpu.fill': 'memory', 'list.bullet.rectangle.portrait.fill': 'fact-check', 'shield.lefthalf.filled': 'security', 'arrow.clockwise': 'refresh', 'account-tree.fill': 'account-tree', 'arrow.triangle.2.circlepath': 'sync', 'waveform.path.ecg': 'timeline', 'checkmark.seal.fill': 'verified', 'doc.text.magnifyingglass': 'find-in-page', 'chevron.right': 'chevron-right', 'lock.fill': 'lock', 'hammer.fill': 'build', 'checkmark': 'check', 'clock.fill': 'schedule', 'person.crop.circle.fill': 'account-circle', 'play.fill': 'play-arrow',
} as const satisfies Record<string, IconName>;

export function IconSymbol({ name, size = 24, color, style, weight }: { name: IconSymbolName; size?: number; color: string | OpaqueColorValue; style?: StyleProp<TextStyle>; weight?: string }) {
  return <MaterialIcons color={color} size={size} name={MAPPING[name]} style={style} />;
}
