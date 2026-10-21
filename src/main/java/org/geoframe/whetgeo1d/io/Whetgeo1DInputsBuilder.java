package org.geoframe.whetgeo1d.io;

import static org.geoframe.whetgeo1d.io.Whetgeo1DInputsHandler.*;

import java.io.File;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;

import org.hortonmachine.dbs.compat.ADb;
import org.hortonmachine.dbs.compat.IHMPreparedStatement;
import org.hortonmachine.dbs.utils.SqlName;
import org.hortonmachine.gears.io.timeseries.OmsTimeSeriesReader;
import org.joda.time.DateTime;
import org.joda.time.format.DateTimeFormat;
import org.joda.time.format.DateTimeFormatter;

/**
 * Writes the WHETGEO-1D input tables read by {@link Whetgeo1DInputsHandler}
 * straight into a database, without passing through CSV files.
 *
 * <p>
 * Collect the pieces of the column with the chainable methods, then
 * {@link #write()} them all at once:
 *
 * <pre>
 * new Whetgeo1DInputsBuilder(db)
 *     .layer(0.0, 50, 0, 1)  // top eta [m], cells, equationStateID, parameterID
 *     .domainBottom(-1.0)
 *     .initialCondition(0.0, 0.0, 285.15)
 *     .initialCondition(-1.0, 0.0, 285.15)
 *     .parameterSet(Map.of("thetaS", 0.38, "thetaR", 0.068, "alpha", 0.8, "n", 1.09))
 *     .vanGenuchtenSlots()
 *     .constantTimeseries(TABLE_TEMPERATURE_TOP_INTERFACE, "2000-01-01 00:00", "2000-03-01 00:00", 60, 20.0)
 *     .write();
 * </pre>
 *
 * <p>
 * The tables are the same that {@link Whetgeo1DInputsHandler#createDbFromCsv(String, ADb)}
 * creates from a folder of CSV files.
 *
 * @author Andrea Antonello (https://g-ant.eu)
 */
public class Whetgeo1DInputsBuilder {

	private static final String LAYER_TYPE = "L";
	private static final int NODATA_INT = (int) NODATA;
	private static final DateTimeFormatter DATE_FORMATTER = DateTimeFormat.forPattern("yyyy-MM-dd HH:mm")
			.withZoneUTC();

	private final ADb db;

	private final List<Object[]> gridRows = new ArrayList<>(); // type, eta, k, equationStateID, parameterID
	private final List<Object[]> pseudoCellRows = new ArrayList<>(); // stackOrder, eqID, parID, cv, note
	private final List<Object[]> icRows = new ArrayList<>(); // eta, psi0, t0
	private final List<Map<String, Double>> parameterSets = new ArrayList<>();
	private final Map<String, String> parameterSlots = new LinkedHashMap<>();
	private final Map<String, Object[]> timeseries = new LinkedHashMap<>(); // table -> {long[], double[]}

	/**
	 * @param db the database to write into, already open. It is not closed by the builder.
	 */
	public Whetgeo1DInputsBuilder(ADb db) {
		this.db = db;
	}

	// ---------------------------------------------------------------- grid

	/**
	 * Add a soil layer, from its top down to the top of the next layer (or to the
	 * {@link #domainBottom(double)}), split into equal cells. Layers can be added in
	 * any order.
	 *
	 * @param etaTop          the elevation of the top of the layer [m], usually 0 at the surface and negative below.
	 * @param cells           the number of cells the layer is split into.
	 * @param equationStateID the index of the closure model of the layer, in the solver's model arrays.
	 * @param parameterID     the parameter set of the layer, 1 for the first {@link #parameterSet(Map)}.
	 * @return the builder.
	 */
	public Whetgeo1DInputsBuilder layer(double etaTop, int cells, int equationStateID, int parameterID) {
		return layer(LAYER_TYPE, etaTop, cells, equationStateID, parameterID);
	}

	/**
	 * Same as {@link #layer(double, int, int, int)}, with an explicit layer type
	 * (stored, not used by the reader).
	 */
	public Whetgeo1DInputsBuilder layer(String type, double etaTop, int cells, int equationStateID,
			int parameterID) {
		if (cells <= 0) {
			throw new IllegalArgumentException("A layer needs at least one cell, got " + cells + " at eta=" + etaTop);
		}
		gridRows.add(new Object[] { type, etaTop, cells, equationStateID, parameterID });
		return this;
	}

	/**
	 * Close the grid: the elevation of the bottom of the lowest layer.
	 *
	 * @param eta the bottom of the domain [m].
	 * @return the builder.
	 */
	public Whetgeo1DInputsBuilder domainBottom(double eta) {
		gridRows.add(new Object[] { LAYER_TYPE, eta, NODATA_INT, NODATA_INT, NODATA_INT });
		return this;
	}

	/**
	 * Add a zero-thickness cell on top of the grid, ex. a surface pond for a
	 * {@code TOP_COUPLED} boundary condition. Cells are stacked in the order they are added.
	 *
	 * @param equationStateID the index of the closure model of the cell, ex. the "Water Depth" one.
	 * @param parameterID     the parameter set of the cell.
	 * @param controlVolume   the storage capacity of the cell [m].
	 * @param note            a free description, or <code>null</code>.
	 * @return the builder.
	 */
	public Whetgeo1DInputsBuilder pseudoCell(int equationStateID, int parameterID, double controlVolume,
			String note) {
		pseudoCellRows.add(new Object[] { pseudoCellRows.size(), equationStateID, parameterID, controlVolume,
				note == null ? "" : note });
		return this;
	}

	// ---------------------------------------------------------------- initial condition

	/**
	 * Add a point of the initial condition. The reader interpolates linearly
	 * between the two topmost points, so give one at the top and one at the bottom
	 * of the domain.
	 *
	 * @param eta  the elevation of the point [m].
	 * @param psi0 the water suction [m], 0 at saturation, negative above the water table.
	 * @param t0   the temperature [K].
	 * @return the builder.
	 */
	public Whetgeo1DInputsBuilder initialCondition(double eta, double psi0, double t0) {
		icRows.add(new Object[] { eta, psi0, t0 });
		return this;
	}

	// ---------------------------------------------------------------- soil parameters

	/**
	 * Add a parameter set: the first one gets parameterID 1, the next 2 and so on.
	 * All sets must have the same names, ex. thetaS, thetaR, n, alpha, Ks, see the
	 * {@code COL_SWRC_*} names of {@link Whetgeo1DInputsHandler}.
	 *
	 * <p>
	 * From groovy, named arguments make it read like a table row:
	 * {@code parameterSet(thetaS: 0.38, thetaR: 0.068, n: 1.09, alpha: 0.8)}.
	 *
	 * <p>
	 * A <code>null</code> value leaves the parameter empty, ex. the {@code psiD} of a
	 * Van Genuchten set in a table that also holds Brooks-Corey sets.
	 *
	 * @param values the values by parameter name.
	 * @return the builder.
	 */
	public Whetgeo1DInputsBuilder parameterSet(Map<String, ?> values) {
		if (values.isEmpty()) {
			throw new IllegalArgumentException("A parameter set needs at least one value.");
		}
		Map<String, Double> set = new LinkedHashMap<>();
		for (Map.Entry<String, ?> entry : values.entrySet()) {
			Object value = entry.getValue();
			if (value != null && !(value instanceof Number)) {
				throw new IllegalArgumentException(
						"Parameter " + entry.getKey() + " must be a number, got: " + value);
			}
			set.put(entry.getKey(), value == null ? null : ((Number) value).doubleValue());
		}
		if (!parameterSets.isEmpty() && !parameterSets.get(0).keySet().equals(set.keySet())) {
			throw new IllegalArgumentException("Parameter set " + (parameterSets.size() + 1) + " has the names "
					+ set.keySet() + ", but the first one has " + parameterSets.get(0).keySet());
		}
		parameterSets.add(set);
		return this;
	}

	/**
	 * Tell which parameter name fills which generic slot of the retention curve
	 * ({@code par1} .. {@code par5}).
	 *
	 * @param name the parameter name, ex. "n".
	 * @param slot the slot, "par1" to "par5".
	 * @return the builder.
	 */
	public Whetgeo1DInputsBuilder parameterSlot(String name, String slot) {
		if (parameterSlots.containsKey(name)) {
			throw new IllegalArgumentException("Parameter name " + name + " is already mapped.");
		}
		parameterSlots.put(name, slot);
		return this;
	}

	/**
	 * The slots of the Van Genuchten retention curve: n is par1, alpha is par2.
	 *
	 * @return the builder.
	 */
	public Whetgeo1DInputsBuilder vanGenuchtenSlots() {
		return parameterSlot(COL_SWRC_N, "par1").parameterSlot(COL_SWRC_ALPHA, "par2");
	}

	/**
	 * The slots of the Brooks-Corey retention curve: n is par1, the air-entry suction
	 * psiD [m] is par2.
	 *
	 * @return the builder.
	 */
	public Whetgeo1DInputsBuilder brooksCoreySlots() {
		return parameterSlot(COL_SWRC_N, "par1").parameterSlot(COL_SWRC_PSI_D, "par2");
	}

	/**
	 * The slots of the Kosugi retention curve: the median suction psiMedian [m] is
	 * par1, sigma is par2.
	 *
	 * @return the builder.
	 */
	public Whetgeo1DInputsBuilder kosugiSlots() {
		return parameterSlot(COL_SWRC_PSI_MEDIAN, "par1").parameterSlot(COL_SWRC_SIGMA, "par2");
	}

	/**
	 * The slots of the Romano bimodal retention curve: the weight w is par1, sigma1
	 * and sigma2 are par2 and par3, the median suctions h1 and h2 [m] are par4 and par5.
	 *
	 * @return the builder.
	 */
	public Whetgeo1DInputsBuilder romanoSlots() {
		return parameterSlot(COL_SWRC_W, "par1").parameterSlot(COL_SWRC_SIGMA1, "par2")
				.parameterSlot(COL_SWRC_SIGMA2, "par3").parameterSlot(COL_SWRC_H1, "par4")
				.parameterSlot(COL_SWRC_H2, "par5");
	}

	// ---------------------------------------------------------------- timeseries

	/**
	 * Add a timeseries, ex. a boundary condition or a meteorological forcing.
	 *
	 * @param table  the table name, ex. {@link Whetgeo1DInputsHandler#TABLE_TEMPERATURE_TOP_INTERFACE}
	 *               or "timeseries_topBC".
	 * @param times  the times, epoch milliseconds UTC.
	 * @param values the values, one per time.
	 * @return the builder.
	 */
	public Whetgeo1DInputsBuilder timeseries(String table, long[] times, double[] values) {
		if (times.length != values.length) {
			throw new IllegalArgumentException(
					"Timeseries " + table + " has " + times.length + " times but " + values.length + " values.");
		}
		if (timeseries.containsKey(table)) {
			throw new IllegalArgumentException("Timeseries " + table + " was already added.");
		}
		timeseries.put(table, new Object[] { times, values });
		return this;
	}

	/**
	 * Add a timeseries with the same value at every step, ex. a fixed boundary temperature.
	 *
	 * @param table       the table name.
	 * @param startDate   the first time, "yyyy-MM-dd HH:mm" UTC.
	 * @param endDate     the last time, included, "yyyy-MM-dd HH:mm" UTC.
	 * @param stepMinutes the time step [min].
	 * @param value       the value.
	 * @return the builder.
	 */
	public Whetgeo1DInputsBuilder constantTimeseries(String table, String startDate, String endDate,
			int stepMinutes, double value) {
		if (stepMinutes <= 0) {
			throw new IllegalArgumentException("The time step must be positive, got " + stepMinutes);
		}
		long start = DATE_FORMATTER.parseDateTime(startDate).getMillis();
		long end = DATE_FORMATTER.parseDateTime(endDate).getMillis();
		if (end < start) {
			throw new IllegalArgumentException("The end date " + endDate + " is before the start date " + startDate);
		}
		long step = stepMinutes * 60_000L;
		int count = (int) ((end - start) / step) + 1;
		long[] times = new long[count];
		double[] values = new double[count];
		for (int i = 0; i < count; i++) {
			times[i] = start + i * step;
			values[i] = value;
		}
		return timeseries(table, times, values);
	}

	/**
	 * Add a timeseries from an OMS timeseries CSV file (the {@code @H} format), which
	 * must have exactly one value column.
	 *
	 * @param table   the table name.
	 * @param csvPath the CSV file.
	 * @return the builder.
	 */
	public Whetgeo1DInputsBuilder timeseriesFromCsv(String table, String csvPath) throws Exception {
		File csv = new File(csvPath);
		String[] headerParts = null;
		for (String line : Files.readAllLines(csv.toPath())) {
			if (line.startsWith("@H,")) {
				headerParts = line.split(",", -1);
				break;
			}
		}
		if (headerParts == null) {
			throw new IllegalArgumentException("No @H header line found in " + csv.getName());
		}
		// headerParts: [@H, timestamp, value_<id>]
		int numValueCols = headerParts.length - 2;
		if (numValueCols != 1) {
			throw new IllegalArgumentException(
					"Expected exactly one value column in " + csv.getName() + ", found " + numValueCols);
		}

		OmsTimeSeriesReader reader = new OmsTimeSeriesReader();
		reader.file = csv.getAbsolutePath();
		reader.read();
		reader.close();

		int size = reader.outData.size();
		long[] times = new long[size];
		double[] values = new double[size];
		int i = 0;
		for (Map.Entry<DateTime, double[]> entry : reader.outData.entrySet()) {
			times[i] = entry.getKey().getMillis();
			values[i] = entry.getValue()[0];
			i++;
		}
		return timeseries(table, times, values);
	}

	// ---------------------------------------------------------------- write

	/**
	 * Check the collected pieces and write them, each into its own table. A table
	 * that does not exist is created, a table that already holds rows is an error.
	 */
	public void write() throws Exception {
		validate();

		if (!parameterSlots.isEmpty()) {
			List<Object[]> rows = new ArrayList<>();
			parameterSlots.forEach((name, parType) -> rows.add(new Object[] { name, parType }));
			insert(TABLE_DICTIONARY, new String[] { COL_DICTIONARY_VAL1 + " TEXT PRIMARY KEY",
					COL_DICTIONARY_VAL2 + " TEXT" }, new String[] { COL_DICTIONARY_VAL1, COL_DICTIONARY_VAL2 },
					rows);
		}
		if (!icRows.isEmpty()) {
			insert(TABLE_IC, new String[] { COL_IC_ETA + " REAL PRIMARY KEY", COL_IC_PSI0 + " REAL",
					COL_IC_T0 + " REAL" }, new String[] { COL_IC_ETA, COL_IC_PSI0, COL_IC_T0 }, icRows);
		}
		if (!gridRows.isEmpty()) {
			// the table is keyed by eta, the column order of the insert is eta first
			List<Object[]> rows = new ArrayList<>();
			for (Object[] r : gridRows) {
				rows.add(new Object[] { r[1], r[0], r[2], r[3], r[4] });
			}
			insert(TABLE_GRID,
					new String[] { COL_GRID_ETA + " REAL PRIMARY KEY", COL_GRID_TYPE + " TEXT",
							COL_GRID_K + " INTEGER", COL_GRID_EQUATION_STATE_ID + " INTEGER",
							COL_GRID_PARAMETER_ID + " INTEGER" },
					new String[] { COL_GRID_ETA, COL_GRID_TYPE, COL_GRID_K, COL_GRID_EQUATION_STATE_ID,
							COL_GRID_PARAMETER_ID },
					rows);
		}
		if (!pseudoCellRows.isEmpty()) {
			insert(TABLE_PSEUDO_CELL,
					new String[] { COL_PSEUDO_CELL_STACK_ORDER + " INTEGER PRIMARY KEY",
							COL_PSEUDO_CELL_EQUATION_STATE_ID + " INTEGER", COL_PSEUDO_CELL_PARAMETER_ID + " INTEGER",
							COL_PSEUDO_CELL_CONTROL_VOLUME + " REAL", COL_PSEUDO_CELL_NOTE + " TEXT" },
					new String[] { COL_PSEUDO_CELL_STACK_ORDER, COL_PSEUDO_CELL_EQUATION_STATE_ID,
							COL_PSEUDO_CELL_PARAMETER_ID, COL_PSEUDO_CELL_CONTROL_VOLUME, COL_PSEUDO_CELL_NOTE },
					pseudoCellRows);
		}
		if (!parameterSets.isEmpty()) {
			List<String> names = new ArrayList<>(parameterSets.get(0).keySet());
			String[] fieldDefs = new String[names.size() + 1];
			String[] columns = new String[names.size() + 1];
			fieldDefs[0] = COL_SWRC_ID + " INTEGER PRIMARY KEY";
			columns[0] = COL_SWRC_ID;
			for (int i = 0; i < names.size(); i++) {
				fieldDefs[i + 1] = names.get(i) + " REAL";
				columns[i + 1] = names.get(i);
			}
			List<Object[]> rows = new ArrayList<>();
			for (int s = 0; s < parameterSets.size(); s++) {
				Object[] row = new Object[columns.length];
				row[0] = s + 1;
				for (int i = 0; i < names.size(); i++) {
					row[i + 1] = parameterSets.get(s).get(names.get(i));
				}
				rows.add(row);
			}
			insert(TABLE_SWRC_PARAMETERS, fieldDefs, columns, rows);
		}
		for (Map.Entry<String, Object[]> entry : timeseries.entrySet()) {
			long[] times = (long[]) entry.getValue()[0];
			double[] values = (double[]) entry.getValue()[1];
			List<Object[]> rows = new ArrayList<>(times.length);
			for (int i = 0; i < times.length; i++) {
				rows.add(new Object[] { times[i], values[i] });
			}
			insert(entry.getKey(), new String[] { COL_TEMPERATURE_TIMESTAMP + " INTEGER PRIMARY KEY",
					COL_TIMESERIES_VALUE + " REAL" }, new String[] { COL_TEMPERATURE_TIMESTAMP, COL_TIMESERIES_VALUE },
					rows);
		}
	}

	private void validate() {
		if (!gridRows.isEmpty()) {
			List<Object[]> sorted = new ArrayList<>(gridRows);
			sorted.sort((a, b) -> Double.compare((double) b[1], (double) a[1]));
			Object[] lowest = sorted.get(sorted.size() - 1);
			if ((int) lowest[2] > 0) {
				throw new IllegalArgumentException(
						"The grid has no bottom: add domainBottom(eta) below the lowest layer (eta=" + lowest[1] + ").");
			}
			if (sorted.stream().noneMatch(r -> (int) r[2] > 0)) {
				throw new IllegalArgumentException("The grid has no layers.");
			}
			LinkedHashSet<Double> etas = new LinkedHashSet<>();
			for (Object[] r : gridRows) {
				if (!etas.add((double) r[1])) {
					throw new IllegalArgumentException("Two grid rows have the same eta=" + r[1]);
				}
			}
			if (!parameterSets.isEmpty()) {
				for (Object[] r : gridRows) {
					int parameterID = (int) r[4];
					if ((int) r[2] > 0 && (parameterID < 1 || parameterID > parameterSets.size())) {
						throw new IllegalArgumentException("The layer at eta=" + r[1] + " uses parameterID "
								+ parameterID + ", but there are " + parameterSets.size() + " parameter sets.");
					}
				}
			}
		}
		if (icRows.size() == 1) {
			throw new IllegalArgumentException(
					"The initial condition needs at least two points, at the top and at the bottom.");
		}
	}

	private void insert(String tableName, String[] fieldDefs, String[] columns, List<Object[]> rows)
			throws Exception {
		SqlName table = SqlName.m(tableName);
		if (!db.hasTable(table)) {
			db.createTable(table, fieldDefs);
		} else if (db.getCount(table) > 0) {
			throw new IllegalStateException("The table " + tableName + " already holds data.");
		}

		String sql = "INSERT INTO " + table.fixedDoubleName + " (" + String.join(", ", columns) + ") VALUES ("
				+ "?, ".repeat(columns.length - 1) + "?)";
		db.execOnConnection(conn -> {
			boolean autoCommit = conn.getAutoCommit();
			conn.setAutoCommit(false);
			try (IHMPreparedStatement pStmt = conn.prepareStatement(sql)) {
				for (Object[] row : rows) {
					for (int i = 0; i < row.length; i++) {
						Object value = row[i];
						if (value == null) {
							pStmt.setObject(i + 1, null);
						} else if (value instanceof Double d) {
							pStmt.setDouble(i + 1, d);
						} else if (value instanceof Integer n) {
							pStmt.setInt(i + 1, n);
						} else if (value instanceof Long n) {
							pStmt.setLong(i + 1, n);
						} else {
							pStmt.setString(i + 1, String.valueOf(value));
						}
					}
					pStmt.addBatch();
				}
				pStmt.executeBatch();
				conn.commit();
			} finally {
				conn.setAutoCommit(autoCommit);
			}
			return null;
		});
	}
}
