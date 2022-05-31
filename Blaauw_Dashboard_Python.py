import panel as pn
pn.extension('tabulator', sizing_mode='stretch_width', notifications=True)
from panel.io.notifications import NotificationArea
NotificationArea.position = 'bottom-right'

import pyvo as vo
import datetime as dt
from astropy.time import Time
import param
from pypika import Table, Criterion, EmptyCriterion, Order
import pandas as pd
from astroquery.simbad import Simbad
from astropy.coordinates import SkyCoord
from astropy import units as u
import warnings
from astroquery.exceptions import TableParseError
from astropy.io import fits
from matplotlib.figure import Figure
from matplotlib import cm
import numpy as np
from pathlib import Path
from astropy.coordinates import Angle

url = "https://vo.astro.rug.nl/tap"
service = vo.dal.TAPService(url)

from pypika.queries import Query, QueryBuilder
from pypika.utils import builder, QueryException
from typing import Any, Optional
from enum import Enum

class Dialects(Enum):
    ADQL = 'adql'

class ADQLQuery(Query):
    @classmethod
    def _builder(cls, **kwargs: Any) -> "QueryBuilder":
        return ADQLQueryBuilder(**kwargs)

class ADQLQueryBuilder(QueryBuilder):
    QUOTE_CHAR = None
    ALIAS_QUOTE_CHAR = '"'
    QUERY_ALIAS_QUOTE_CHAR = ''
    QUERY_CLS = ADQLQuery

    def __init__(self, **kwargs: Any) -> None:
        super().__init__(dialect=Dialects.ADQL, **kwargs)
        self._top = None
        self._contains_point = None
        self._contains_box = None

    @builder
    def top(self, value):
        try:
            self._top = int(value)
        except ValueError:
            raise QueryException('TOP value must be an integer')

    def get_sql(self, *args, **kwargs):
        return super(ADQLQueryBuilder, self).get_sql(*args, groupby_alias=False, **kwargs)

    def _top_sql(self):
        if self._top:
            return 'TOP {} '.format(self._top)
        else:
            return ''

    @builder
    def contains_point(self, ra, dec, from_table, cone_radius):
        self._contains_point = [ra, dec, from_table, cone_radius]

    def _contains_point_sql(self):
        if self._contains_point:
            return "CONTAINS(POINT('ICRS', {}, {}), CIRCLE('ICRS', {}.ra, {}.dec, {}))=1 AND ".format(
                self._contains_point[0], self._contains_point[1], self._contains_point[2], self._contains_point[2], self._contains_point[3])
        else:
            return ''

    @builder
    def contains_box(self, ra, dec, from_table, width, height):
        self._contains_box = [ra, dec, from_table, width, height]

    def _contains_box_sql(self):
        if self._contains_box:
            return "CONTAINS(POINT('ICRS', {}, {}), BOX('ICRS', {}.ra, {}.dec, {}, {}))=1 AND ".format(
                self._contains_box[0], self._contains_box[1], self._contains_box[2], self._contains_box[2], self._contains_box[3], self._contains_box[4])
        else:
            return ''

    def _select_sql(self, **kwargs: Any) -> str:
        return "SELECT {distinct}{top}{select}".format(
            top=self._top_sql(),
            distinct="DISTINCT " if self._distinct else "",
            select=",".join(term.get_sql(with_alias=True, subquery=True, **kwargs) for term in self._selects)
        )

    @builder
    def _where_sql(self, quote_char: Optional[str] = None, **kwargs: Any) -> str:
        return " WHERE {contains_point}{contains_box}{where}".format(
            contains_point=self._contains_point_sql(),
            contains_box = self._contains_box_sql(),
            where=self._wheres.get_sql(quote_char=quote_char, subquery=True, **kwargs)
        )

home_dir = '../Data/Raw/'

options = {'Filename': 'filename',
           'Observation date': 'kw_DATE_OBS',
           'Image type': 'kw_IMAGETYP',
           'Right ascension': 'ra',
           'Declination': 'dec',
           'Object': 'kw_OBJECT',
           'Filter': 'kw_FILTER'}
options_reverse = {'filename': 'Filename',
                   'kw_DATE_OBS': 'Observation date',
                   'kw_IMAGETYP': 'Image type',
                   'ra': 'Right ascension',
                   'dec': 'Declination',
                   'kw_OBJECT': 'Object',
                   'kw_FILTER': 'Filter'}
class CompileQuery(param.Parameterized):
    begin_date = dt.datetime(2008, 5, 1, 12)
    end_date = dt.datetime(2020, 4, 23, 12)
    dates = param.DateRange((begin_date, end_date), precedence=1)
    dates_order = param.Selector({'Descending': Order.desc, 'Ascending': Order.asc}, precedence=1)
    reset_date_button = param.Action(lambda x: x.param.trigger('reset_date_button'), label='Reset date', precedence=1)
    select_columns = param.ListSelector(default=['Filename',
                                                 'Image type',
                                                 'Filter',
                                                 'Observation date',
                                                 'Right ascension',
                                                 'Declination',
                                                 'Object'],
                                        objects=options.keys(),
                                        precedence=1)
    ra_dec_not_null = param.Boolean(default=True, precedence=1)
    object_or_coordinates = param.Selector(['Object', 'Coordinates'], precedence=1)
    select_object = param.String(default='', precedence=1)
    search_coordinates = param.String(default='' , precedence=-1)
    box_or_cone = param.Selector(['Box', 'Cone'], precedence=1)
    cone_radius = param.Number(default=20, bounds=(5, 40), precedence=-1)
    box_width = param.Number(1700, precedence=1)
    box_height = param.Number(1100, precedence=1)
    nr_entries = param.Number(100, precedence=1)

    @param.depends('box_or_cone', watch=True)
    def box_cone_search(self):
        if self.box_or_cone == 'Box':
            self.param.cone_radius.precedence = -1
            self.param.box_width.precedence = 1
            self.param.box_height.precedence = 1
        else:
            self.param.cone_radius.precedence = 1
            self.param.box_width.precedence = -1
            self.param.box_height.precedence = -1

    @param.depends('object_or_coordinates', watch=True)
    def object_coordinates_search(self):
        if self.object_or_coordinates == 'Object':
            self.param.select_object.precedence = 1
            self.param.search_coordinates.precedence = -1
        else:
            self.select_object = ''
            self.param.select_object.precedence = -1
            self.param.search_coordinates.precedence = 1

    @param.depends('reset_date_button', watch=True)
    def reset_date(self):
        self.dates = (dt.datetime(2008, 5, 1, 12), dt.datetime(2022, 5, 1, 12))

    def table_string(self):
        table_name = 'observations.raw'
        return table_name

    def from_query(self):
        raw_observations = Table(self.table_string())
        return raw_observations

    @param.depends('select_columns')
    def select_query(self):
        select_options = ''
        from_table = self.from_query()
        for i in self.select_columns:
            select_options += str(from_table).strip('"') + '.' + options[i] + ','
        return select_options.rstrip(',')

    @param.depends('select_object')
    def find_object(self):
        if self.select_object == '':
            return True
        else:
            object_name = self.select_object
            with warnings.catch_warnings():
                warnings.filterwarnings('error')
                try:
                    simbad_object = Simbad.query_object(object_name)
                    coordinates = SkyCoord(simbad_object['RA'], simbad_object['DEC'], unit=(u.hourangle, u.deg))
                    object_ra, object_dec = coordinates.ra.deg[0], coordinates.dec.deg[0]
                    return object_ra, object_dec
                except (Warning, TableParseError):
                    return False

    @param.depends('search_coordinates')
    def coordinates_to_degrees(self):
        if self.search_coordinates == '':
            return
        right_ascension, declination = self.search_coordinates.split(',')
        if len(right_ascension) > 1:
            if float(right_ascension.split(' ')[0]) >= 24:
                right_ascension = Angle(right_ascension, unit=u.degree)
                right_ascension.wrap_at(360 * u.degree, inplace=True)
            else:
                right_ascension = Angle(right_ascension, unit=u.hourangle)
                right_ascension.wrap_at(24 * u.hourangle, inplace=True)
        else:
            if float(right_ascension) >= 24:
                right_ascension = Angle(right_ascension, unit=u.degree)
                right_ascension.wrap_at(360 * u.degree, inplace=True)
            else:
                right_ascension = Angle(right_ascension, unit=u.hourangle)
                right_ascension.wrap_at(24 * u.hourangle, inplace=True)
        declination = Angle(declination, unit=u.deg)
        right_ascension = right_ascension.deg
        declination = declination.deg
        if declination < -90:
            while declination < -90:
                declination += 90
        elif declination > 90:
            while declination > 90:
                declination -= 90
        return right_ascension, declination

    @param.depends('ra_dec_not_null')
    def not_null_check(self):
        if self.ra_dec_not_null:
            return True

    @param.depends('dates')
    def date_query(self):
        start_date = Time(self.dates[0]).jd
        end_date = Time(self.dates[1]).jd
        return start_date, end_date

class FetchData(CompileQuery):
    search_object_button = param.Action(lambda x: x.param.trigger('search_object_button'), label='Search', precedence=5)
    reset_search_button = param.Action(lambda x: x.param.trigger('reset_search_button'), label='Reset search', precedence=5)

    @param.depends('reset_search_button', watch=True)
    def reset_search(self):
        self.dates = (self.begin_date, self.end_date)
        self.ra_dec_not_null = True
        self.select_object = ''
        self.cone_radius = 20
        self.param.trigger('search_object_button')

    @param.depends('search_object_button')
    def data(self):
        if not self.find_object():
            return
        observation_table = self.from_query()
        dates = self.date_query()
        query = ADQLQuery.from_(
            observation_table
        ).select(
            self.select_query()
        ).where(
            Criterion.all([
                observation_table.ra.isnotnull() if self.not_null_check() else EmptyCriterion(),
                observation_table.dec.isnotnull() if self.not_null_check() else EmptyCriterion(),
                observation_table.obs_jd[dates[0]:dates[1]]
            ])
        ).orderby(observation_table.obs_jd, order=self.dates_order).top(self.nr_entries)
        if self.select_object != '':
            right_ascension, declination = self.find_object()
            if self.box_or_cone == 'Cone':
                query = query.contains_point(right_ascension, declination, self.table_string(), self.cone_radius/60)
            elif self.box_or_cone == 'Box':
                query = query.contains_box(self.find_object()[0], self.find_object()[1], self.table_string(), self.box_width/3600, self.box_height/3600)
        elif self.select_object == '' and self.search_coordinates != '':
            right_ascension, declination = self.coordinates_to_degrees()
            if self.box_or_cone == 'Cone':
                query = query.contains_point(right_ascension, declination, self.table_string(), self.cone_radius/60)
            elif self.box_or_cone == 'Box':
                query = query.contains_box(right_ascension, declination, self.table_string(), self.box_width/3600, self.box_height/3600)
        result = service.search(query.get_sql())
        data_pandas = result.to_table().to_pandas()
        data_pandas['ra'] = data_pandas['ra'].apply(self.ra_to_hms)
        data_pandas['dec'] = data_pandas['dec'].apply(self.dec_to_dms)
        return data_pandas

    def ra_to_hms(self, nr_angle):
        angle = Angle(nr_angle, u.degree)
        return angle.to_string(unit=u.hour, sep=('h:', 'm:', 's'), precision=3)

    def dec_to_dms(self, nr_angle):
        angle = Angle(nr_angle, u.degree)
        return angle.to_string(unit=u.degree, sep=('d:', 'm:', 's'), precision=3)

    @param.depends('search_object_button')
    def display_object_error(self):
        if not self.find_object():
            message = "SIMBAD could not find the object " + str(self.select_object) + " that you searched for."
            pn.state.notifications.error(message, duration=0)

    @param.depends('search_object_button', 'data')
    def display_data_error(self):
        if self.find_object() and self.data().values.size == 0:
            message = "No entries could be found with the options you specified."
            pn.state.notifications.error(message, duration=0)

class CreateTable(FetchData):
    table = param.DataFrame(columns=CompileQuery().select_columns, precedence=-1)
    if FetchData().data() is None:
        table.default = pd.DataFrame(columns=CompileQuery().select_columns)
    else:
        table.default = FetchData().data()
        table.columns = options_reverse

    def __init__(self, **params):
        super().__init__(**params)
        self.table_widget = pn.widgets.Tabulator(self.table,
                                                 disabled=True,
                                                 pagination='local',
                                                 show_index=False,
                                                 layout='fit_columns',
                                                 selectable='checkbox-single',
                                                 text_align='left',
                                                 min_height=500,
                                                 titles=options_reverse,
                                                 selection=[8, 14])

    @param.depends('data', watch=True)
    def full_table(self):
        if self.data() is None:
            self.table_widget.value = pd.DataFrame(columns=self.select_columns)
        else:
            self.table_widget.value = self.data()
            self.table_widget.titles = options_reverse

    selected_table_widget = pn.widgets.Tabulator(pd.DataFrame(),
                                                 disabled=True,
                                                 pagination='local',
                                                 page_size=11,
                                                 show_index=False,
                                                 layout='fit_data',
                                                 selectable='checkbox',
                                                 text_align='left',
                                                 min_height=350,
                                                 titles=options_reverse)

    @param.depends('table_widget.selection')
    def selected_table(self):
        selection = self.table_widget.selected_dataframe
        if selection.values.size == 0:
            self.selected_table_widget.value = pd.DataFrame(columns=['filename', 'kw_FILTER', 'kw_OBJECT'])
        else:
            selection = selection[['filename', 'kw_FILTER', 'kw_OBJECT']]
            selection['filename'] = selection['filename'].apply(self.give_file_name)
            self.selected_table_widget.value = selection
        return self.selected_table_widget

    def give_file_name(self, filepath):
        return Path(filepath).name


class CreatePlots(CreateTable):
    plot_button = param.Action(lambda x: x.param.trigger('plot_button'), label='Plot', precedence=5)

    def __init__(self, **params):
        super().__init__(**params)
        self.cannot_find_file = False

    def get_file_list(self):
        selected = self.selected_table_widget.selected_dataframe
        if len(selected) == 0 or len(selected) > 3:
            return None
        filename_list = []
        for i in selected['filename']:
            filepath = home_dir + i
            if Path(filepath).is_file():
                filename_list.append(home_dir + i)
            else:
                self.cannot_find_file = True
                return None
        return filename_list

    @param.depends('plot_button')
    def single_plot(self, filename):
        fits_file = fits.open(filename)
        fits_data = fits_file[0].data
        filter_type = fits_file[0].header['FILTER']
        fits_file.close()
        vmin = np.percentile(fits_data, 5)
        vmax = np.percentile(fits_data, 95)
        fig = Figure(figsize=(9, 6))
        ax = fig.subplots()
        short_name = self.give_file_name(filename)
        ax.set_title(short_name + " (filter: " + filter_type + ")")
        img = ax.imshow(fits_data, interpolation='none', origin='lower', cmap=cm.gray, vmin=vmin, vmax=vmax)
        fig.colorbar(img, shrink=0.8)
        img.axes.get_xaxis().set_visible(False)
        img.axes.get_yaxis().set_visible(False)
        return pn.pane.Matplotlib(fig, tight=True, sizing_mode='scale_both')

    @param.depends('plot_button')
    def file_not_found(self):
        if self.cannot_find_file:
            self.cannot_find_file = False
            message = "One of the files you tried to plot could not be found."
            pn.state.notifications.error(message, duration=0)

    @param.depends('plot_button')
    def plot1(self):
        if self.get_file_list() is None:
            return pn.Card()
        else:
            file = self.get_file_list()[0]
            return pn.Card(self.single_plot(file))

    @param.depends('plot_button')
    def plot2(self):
        if self.get_file_list() is None:
            return pn.Card()
        elif len(self.get_file_list()) < 2:
            return pn.Card()
        else:
            file = self.get_file_list()[1]
            return pn.Card(self.single_plot(file))

    @param.depends('plot_button')
    def plot3(self):
        if self.get_file_list() is None:
            return pn.Card()
        elif len(self.get_file_list()) < 3:
            return pn.Card()
        else:
            file = self.get_file_list()[2]
            return pn.Card(self.single_plot(file))

    @param.depends('plot_button')
    def selection_error(self):
        if len(self.selected_table_widget.selection) == 0:
            message = "You have to select at least one file to plot."
            pn.state.notifications.error(message, duration=0)
        elif len(self.selected_table_widget.selection) > 3:
            message = "You can only select three files at the same time for plotting."
            pn.state.notifications.error(message, duration=0)


class Statistics(CreatePlots):
    # years = [2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2020]
    years = ["2009", "2010", "2011", "2012", "2013", "2014", "2015", "2016", "2017", "2018", "2019", "2020"]
    night = param.Date(dt.date(2020, 4, 22), precedence=1)
    specific_year = param.ListSelector(default=years, objects=years, precedence=1)
    statistics_button = param.Action(lambda x: x.param.trigger('statistics_button'), label='Fetch statistics', precedence=5)

    def __init__(self, **params):
        super().__init__(**params)
        self.all_data = pd.DataFrame()
        self.bias_df = pd.DataFrame()
        self.dark_df = pd.DataFrame()
        self.flat_df = pd.DataFrame()

    def gather_data(self):
        begin_night = dt.datetime.combine(self.night, dt.time(12, 0))
        begin_night_jd = Time(begin_night).jd
        end_night = begin_night + dt.timedelta(days=1)
        end_night_jd = Time(end_night).jd
        table_name = self.from_query()
        query = ADQLQuery.from_(
            table_name
        ).select(
            table_name.filename, table_name.kw_IMAGETYP, table_name.kw_XBINNING
        ).where(
            table_name.obs_jd[begin_night_jd:end_night_jd]
        ).top(3000)
        result = service.search(query.get_sql())
        self.all_data = result.to_table().to_pandas()

    @param.depends('statistics_button', watch=True)
    def calibration_files(self):
        self.gather_data()
        self.bias_df = self.all_data[self.all_data['kw_IMAGETYP']=='Bias Frame']
        self.dark_df = self.all_data[self.all_data['kw_IMAGETYP']=='Dark Frame']
        self.flat_df = self.all_data[self.all_data['kw_IMAGETYP']=='Flat Field']

    def bias_statistics(self, binning):
        bias_median = np.array([])
        read_noise = np.array([])
        bias_bin_df = self.bias_df[self.bias_df['kw_XBINNING']==binning]
        if bias_bin_df.empty:
            return []
        else:
            for i in bias_bin_df['filename']:
                filename = self.give_file_name(i)
                filepath = home_dir + filename
                fits_file = fits.open(filepath)
                fits_data = fits_file[0].data
                fits_file.close()
                bias_median = np.append(bias_median, np.median(fits_data))
                read_noise = np.append(read_noise, np.var(fits_data)**0.5)
            bias_merged = np.column_stack((bias_median, read_noise))
            return bias_merged

    def dark_statistics(self, binning):
        dark_median = np.array([])
        dark_exposure = np.array([])
        dark_bin_df = self.dark_df[self.dark_df['kw_XBINNING']==binning]
        if dark_bin_df.empty:
            return []
        else:
            for i in self.dark_df['filename']:
                filename = self.give_file_name(i)
                filepath = home_dir + filename
                fits_file = fits.open(filepath)
                fits_data = fits_file[0].data
                exposure_time = fits_file[0].header['EXPTIME']
                fits_file.close()
                dark_median = np.append(dark_median, np.median(fits_data))
                dark_exposure = np.append(dark_exposure, exposure_time)
            dark_merged = np.column_stack((dark_median, dark_exposure))
            return dark_merged

    def flat_statistics(self, binning):
        flat_median = np.array([])
        flat_filters = np.array([])
        flat_bin_df = self.flat_df[self.flat_df['kw_XBINNING']==binning]
        if flat_bin_df.empty:
            return []
        else:
            for i in self.flat_df['filename']:
                filename = self.give_file_name(i)
                filepath = home_dir + filename
                fits_file = fits.open(filepath)
                fits_data = fits_file[0].data
                file_filter = fits_file[0].header['FILTER']
                fits_file.close()
                flat_median = np.append(flat_median, np.median(fits_data))
                flat_filters = np.append(flat_filters, file_filter)
            flat_merged = np.column_stack((flat_median, flat_filters))
            return flat_merged

colours = ['g', 'r', 'm', 'c']
linestyles = ['dotted', 'dashed', 'dashdot', 'solid']

class StatisticsPlots(Statistics):
    @param.depends('statistics_button')
    def bias_median_plot(self):
        if self.bias_df.empty:
            return pn.Card()
        else:
            fig = Figure(figsize=(9, 6))
            ax = fig.subplots()
            ax.set_title("Median of the bias frames")
            binning_list = self.bias_df.kw_XBINNING.unique()
            binning_list = np.sort(binning_list)
            for i in binning_list:
                bias_median = self.bias_statistics(binning=i)[:, 0]
                binning_string = " (" + str(i) + "x" + str(i) + "binning)"
                if bias_median.size != 0:
                    ax.axvline(np.median(bias_median), label="Median of the night" + binning_string, linewidth=2, linestyle=linestyles[i-1], color=colours[i-1], alpha=0.5)
                old_bias_median = np.loadtxt('merged_statistics/bias_bin' + str(i) + '/bias_bin' + str(i) + '.txt', dtype='str')[:, 1].astype(float)
                vmin = np.percentile(old_bias_median, 10)
                vmax = np.percentile(old_bias_median, 90)
                ax.hist(old_bias_median, bins=int(len(old_bias_median)/5), label="Previous data" + binning_string, alpha=0.75, range=(vmin, vmax))
            ax.legend()
            return pn.Card(pn.pane.Matplotlib(fig, tight=True, sizing_mode='scale_both'))

    @param.depends('statistics_button')
    def read_noise_plot(self):
        if self.bias_df.empty:
            return pn.Card()
        else:
            fig = Figure(figsize=(9, 6))
            ax = fig.subplots()
            ax.set_title("Median of the read noise")
            binning_list = self.bias_df.kw_XBINNING.unique()
            binning_list = np.sort(binning_list)
            for i in binning_list:
                read_noise = self.bias_statistics(binning=i)[:, 1]
                binning_string = " (" + str(i) + "x" + str(i) + "binning)"
                if read_noise.size != 0:
                    ax.axvline(np.median(read_noise), label="Median of the night" + binning_string, linewidth=2, linestyle=linestyles[i-1], color=colours[i-1], alpha=0.5)
                old_read_noise = np.loadtxt('merged_statistics/bias_bin' + str(i) + '/bias_bin' + str(i) + '.txt', dtype='str')[:, 2].astype(float)
                vmin = np.percentile(old_read_noise, 10)
                vmax = np.percentile(old_read_noise, 90)
                ax.hist(old_read_noise, bins=int(len(old_read_noise)/5), label="Previous data" + binning_string, alpha=0.75, range=(vmin, vmax))
            ax.legend()
            return pn.Card(pn.pane.Matplotlib(fig, tight=True, sizing_mode='scale_both'))

    @param.depends('statistics_button')
    def dark_median_plot(self):
        if self.dark_df.empty:
            return pn.Card()
        else:
            fig = Figure(figsize=(9, 6))
            ax = fig.subplots()
            ax.set_title("Median of the dark frames")
            binning_list = self.dark_df.kw_XBINNING.unique()
            binning_list = np.sort(binning_list)
            for i in binning_list:
                dark_median = self.dark_statistics(binning=i)[:, 0]
                binning_string = " (" + str(i) + "x" + str(i) + "binning)"
                if dark_median.size != 0:
                    ax.axvline(np.median(dark_median), label="Median of the night" + binning_string, linewidth=2, linestyle=linestyles[i-1], color=colours[i-1], alpha=0.5)
                old_dark_median = np.loadtxt('merged_statistics/dark_bin' + str(i) + '/dark_bin' + str(i) + '.txt', dtype='str')[:, 1].astype(float)
                vmin = np.percentile(old_dark_median, 10)
                vmax = np.percentile(old_dark_median, 90)
                ax.hist(old_dark_median, bins=int(len(old_dark_median)/5), label="Previous data" + binning_string, alpha=0.75, range=(vmin, vmax))
            ax.legend()
            return pn.Card(pn.pane.Matplotlib(fig, tight=True, sizing_mode='scale_both'))

    @param.depends('statistics_button')
    def flat_median_plot(self):
        if self.flat_df.empty:
            return pn.Card()
        else:
            fig = Figure(figsize=(9, 6))
            ax = fig.subplots()
            ax.set_title("Median of the flat frames")
            binning_list = self.flat_df.kw_XBINNING.unique()
            binning_list = np.sort(binning_list)
            for i in binning_list:
                flat_median = self.flat_statistics(binning=i)[:, 0].astype(float)
                binning_string = " (" + str(i) + "x" + str(i) + "binning)"
                if flat_median.size != 0:
                    ax.axvline(np.median(flat_median), label="Median of the night" + binning_string, linewidth=2, linestyle=linestyles[i-1], color=colours[i-1], alpha=0.5)
                old_flat_median = np.loadtxt('merged_statistics/flat_bin' + str(i) + '/flat_bin' + str(i) + '.txt', dtype='str')[:, 1].astype(float)
                vmin = np.percentile(old_flat_median, 10)
                vmax = np.percentile(old_flat_median, 90)
                ax.hist(old_flat_median, bins=int(len(old_flat_median)/5), label="Previous data" + binning_string, alpha=0.75, range=(vmin, vmax))
            ax.legend()
            return pn.Card(pn.pane.Matplotlib(fig, tight=True, sizing_mode='scale_both'))


class CustomGrid(pn.GridBox):
    def __init__(self, *objects, **params):
        super().__init__(*objects, **params, ncols=2, nrows=1)

class CreateView(StatisticsPlots):
    def panel(self):
        widgets_primary = {
            'dates': {'widget_type': pn.widgets.DatetimeRangePicker, 'max_width': 350},
            'dates_order': {'widget_type': pn.widgets.Select, 'max_width': 160},
            'reset_date_button': {'widget_type': pn.widgets.Button, 'button_type': 'warning', 'max_width': 160, 'align': 'end'},
            'select_columns': {'widget_type': pn.widgets.CrossSelector, 'definition_order': False},
            'ra_dec_not_null': pn.widgets.Checkbox,
            'object_or_coordinates': {'widget_type': pn.widgets.RadioBoxGroup, 'inline': True},
            'select_object': pn.widgets.TextInput,
            'search_coordinates': {'widget_type': pn.widgets.TextInput, 'name': 'Coordinates (hh mm ss.ms, dd mm ss.ms)'},
            'box_or_cone': {'widget_type': pn.widgets.RadioBoxGroup, 'inline': True},
            'cone_radius': {'widget_type': pn.widgets.FloatSlider, 'name': 'Cone radius (arc minutes)' , 'step': 0.5},
            'box_width': {'widget_type': pn.widgets.FloatInput, 'name': 'Box width (arc seconds)', 'max_width': 160},
            'box_height': {'widget_type': pn.widgets.FloatInput, 'name': 'Box height (arc seconds)', 'max_width': 160},
            'nr_entries': {'widget_type': pn.widgets.FloatInput, 'name': '# of entries'},
            'search_object_button': {'widget_type': pn.widgets.Button, 'button_type': 'primary'},
            'reset_search_button': {'widget_type': pn.widgets.Button, 'button_type': 'warning'}
        }
        settings_primary1 = pn.Row(
            pn.Param(self, widgets=widgets_primary, width=385, sizing_mode="fixed", name="Settings", parameters=[
                'dates'
            ])
        )
        settings_primary2 = pn.Row(
            pn.Param(self, widgets=widgets_primary, width=385, sizing_mode="fixed", default_layout=CustomGrid, show_name=False, parameters=[
                'dates_order',
                'reset_date_button'
            ])
        )
        settings_primary3 = pn.Row(
            pn.Param(self, widgets=widgets_primary, width=385, sizing_mode="fixed", show_name=False, parameters=[
                'select_columns',
                'object_or_coordinates',
                'select_object',
                'search_coordinates',
                'box_or_cone',
                'cone_radius'
            ])
        )
        settings_primary4 = pn.Row(
            pn.Param(self, widgets=widgets_primary, width=385, sizing_mode="fixed", default_layout=CustomGrid, show_name=False, parameters=[
                'box_width',
                'box_height'
            ])
        )
        settings_primary5 = pn.Row(
            pn.Param(self, widgets=widgets_primary, width=385, sizing_mode="fixed", show_name=False, parameters=[
                'nr_entries',
                'search_object_button',
                'reset_search_button'
            ])
        )
        settings_tabs = pn.Tabs(
            ('Query', pn.Column(
                settings_primary1,
                settings_primary2,
                settings_primary3,
                settings_primary4,
                settings_primary5
            ))
        )

        widgets_secondary = {
            'night': {'widget_type': pn.widgets.DatePicker},
            'specific_year' : {'widget_type': pn.widgets.MultiSelect, 'height': 215},
            'statistics_button': {'widget_type': pn.widgets.Button, 'button_type': 'primary'}
        }
        settings_secondary = pn.Row(
            pn.Param(self, widgets=widgets_secondary, width=385, sizing_mode="fixed", name="Settings", parameters=[
                'night',
                'specific_year',
                'statistics_button'
            ])
        )
        settings_tabs.append(
            ('Statistics', pn.Column(
                settings_secondary
            ))
        )

        bootstrap.sidebar.append(settings_tabs)

        panel_plot_button = {'plot_button': {'widget_type': pn.widgets.Button, 'button_type': 'primary', 'name': 'Plot'}}
        plot_grid = pn.Column(
            pn.Row(
                pn.Card(self.selected_table,
                        pn.Param(self, widgets=panel_plot_button, show_name=False, parameters=['plot_button'])),
                self.plot1
            ),
            pn.Row(
                self.plot2,
                self.plot3
            ),
            self.selection_error,
            self.file_not_found
        )

        statistics_grid = pn.Column(
            pn.Row(
                self.bias_median_plot,
                self.read_noise_plot
            ),
            pn.Row(
                self.dark_median_plot,
                self.flat_median_plot
            )
        )

        main_tabs = pn.Tabs(
            ('Data Table', pn.Column(self.table_widget,
                                     self.display_object_error,
                                     self.display_data_error)
             ),
            ('Data Plotting', plot_grid),
            ('Statistics Plotting', statistics_grid)
        )

        def change_settings_tabs(target, event):
            if event.new == 0 or event.new == 1:
                settings_tabs.active = 0
            else:
                settings_tabs.active = 1

        main_tabs.link(settings_tabs, callbacks={'active': change_settings_tabs})

        bootstrap.main.append(main_tabs)
        return bootstrap


bootstrap = pn.template.BootstrapTemplate(title='Bootstrap Template', sidebar_width=400)

test = CreateView()
test.panel().servable()