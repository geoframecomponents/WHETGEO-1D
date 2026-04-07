package org.geoframe.whetgeo1d.utils;

/**
 * Helper class for enums.
 * 
 * @author Andrea Antonello
 */
public final class EnumUtils {

    private EnumUtils() {
    }

    public static <E extends Enum<E>> E fromString(Class<E> enumClass, String value) {
        if (value == null) {
            throw new IllegalArgumentException("Value cannot be null");
        }

        String normalized = normalize(value);

        for (E constant : enumClass.getEnumConstants()) {
            if (normalize(constant.name()).equals(normalized)
                    || normalize(constant.toString()).equals(normalized)) {
                return constant;
            }
        }

        throw new IllegalArgumentException("Unknown " + enumClass.getSimpleName() + ": " + value);
    }

    private static String normalize(String value) {
        return value.trim()
                .replace("_", "")
                .replace("-", "")
                .replace(" ", "")
                .toLowerCase();
    }
}